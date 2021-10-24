//SuperTrack
#include "testCUDA.cuh"
#include "utils.hh"
//ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TEntryList.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TMath.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TPad.h"
#include "ROOT/TProcessExecutor.hxx"
//STD
#include <unordered_map>
#include <vector>
#include <iterator>
#include <tuple>
#include <filesystem>
//CUDA Libraries
#include <cuda.h>
#include <curand.h>
//Thrust
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
//CUB (Cuda UnBound)
#include <cub/cub.cuh>

//
//Structs and typedefs
//

struct Track
{ 
	double x;
	double y;
	double z;
	double edep; 
};
typedef struct Track Track; //this is a typdef which maps from struct Track --> Track. Just saves us from writing struct when we refer to the struct later.

struct SphericalGeometry
{
	double greatestSphereOffset;
	double sphereRadius;
	long numSpheresLinear;
};
typedef struct SphericalGeometry SphericalGeometry;



void LoadTrack(const std::tuple<Int_t,Int_t,Int_t,TString> &input, Track **hostTrack, Track **deviceTrack)
{
	//Open the file in each process and make a Tree Reader
	TFile f = TFile(std::get<3>(input));
	TTreeReader trackReader("Tracks", &f);
	trackReader.SetEntriesRange(std::get<0>(input),std::get<1>(input));
	TTreeReaderValue<double_t> xReader(trackReader, "x [nm]");
	TTreeReaderValue<double_t> yReader(trackReader, "y [nm]");
	TTreeReaderValue<double_t> zReader(trackReader, "z [nm]");
	TTreeReaderValue<double_t> edepReader(trackReader, "edep [eV]");

	std::cout << "thread #: " << std::get<2>(input) << " starting at: " << std::to_string(std::get<0>(input)) << std::endl;

	//Determine size of arrays. Define them. Then allcate unified memory on CPU and GPU
	long nVals = std::get<1>(input) - std::get<0>(input) + 1; //+1 because number of values includes first and last value
	size_t trackSize = nVals * sizeof(double);
	size_t trackStructSize = nVals*sizeof(Track);

	//malloc and cudaMalloc our arrays respectively
	*hostTrack = (Track *)malloc(trackStructSize);
	cudaMalloc(deviceTrack,trackStructSize);

	//Fill the unified memory arrays from the CPU
	for (long loopnum = 0; trackReader.Next(); loopnum++) 
	{
		(*hostTrack)[loopnum].x = *xReader;
		(*hostTrack)[loopnum].y = *yReader;
		(*hostTrack)[loopnum].z = *zReader;
		(*hostTrack)[loopnum].edep = *edepReader;
	}

	//Copy track to GPU memory
	cudaMemcpy(*deviceTrack,*hostTrack,trackStructSize,cudaMemcpyHostToDevice);

	//TODO: free host track, maybe we can't free it yet because it's still being copied right
}

void GenerateRandomXYShift(const std::tuple<Int_t,Int_t,Int_t,TString> &input, float **randomVals, const int &nSamples, const long &random_seed)
{
	cudaMalloc(randomVals,2*sizeof(float)*nSamples); 
	
	//Random number generation on GPU
	curandGenerator_t randGenerator;
	curandCreateGenerator(&randGenerator,CURAND_RNG_PSEUDO_DEFAULT); //consider changing this to Mersenne Twister later
	curandSetPseudoRandomGeneratorSeed(randGenerator,random_seed+std::get<2>(input));
	curandGenerateUniform(randGenerator,*randomVals,2*nSamples);
	curandDestroyGenerator(randGenerator);
	cudaDeviceSynchronize();
}

void GenerateLogHistogram(double **logBins, int **histogramVals, int **histogramValsAccumulated, int nbins, float binLowerMagnitude, float binUpperMagnitude)
{
	//Get the device Id for active GPU
	int deviceId;
	cudaGetDevice(&deviceId);   

	//Fill the log bins and send to the device
	cudaMallocManaged(logBins, (nbins+1)*sizeof(double));
	LogSpace(binLowerMagnitude,binUpperMagnitude,nbins,*logBins);
	cudaMemPrefetchAsync(*logBins,(nbins+1)*sizeof(double),deviceId);

	//TODO: Change to unmanaged memory later
	cudaMallocManaged(histogramVals,nbins*sizeof(int));
	cudaMallocManaged(histogramValsAccumulated,nbins*sizeof(int));

	//Set arrays to zero
	cudaMemset(*histogramVals,0,nbins*sizeof(int));
	cudaMemset(*histogramValsAccumulated,0,nbins*sizeof(int));
}

__global__ void FilterInScoringBox(double greatestSphereOffset, double sphereRadius, long numSpheresLinear, float* randomVals, Track *inputTrack, Track *outputTrack, int numElements, int *numElementsCompacted, int oversampleIterationNumber)
{
	//This function
	//1.) Applies the random shift to the x,y coordinates
	//2.) Checks which edep events are in the box
	//3.) Performs stream compaction on those events which are in the box
	//We are using stream compaction to avoid a monolithic kernel with large blocks within if-statements which reduces warp efficiency

	//Determine index and stride
 	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;
	//printf("index %d, stride %d \n",index,stride);
	
	//Convert random shifts in to appropriate range
	double x_shift = ((randomVals[(oversampleIterationNumber*2)]*greatestSphereOffset*2)-greatestSphereOffset);
	double y_shift = ((randomVals[(oversampleIterationNumber*2+1)]*greatestSphereOffset*2)-greatestSphereOffset);
	//Put definitions outside of for-loop to prevent repeat constructor calls
	double x_shifted; double y_shifted; int outputIndex;
	//The value we compare to, to check if it's in the box
	double box_edge = abs(greatestSphereOffset)+(sphereRadius);

	//Loop over all the energy deposition points
	for (int i = index; i < numElements; i+=stride)
	{
		//Apply random shift
		x_shifted = inputTrack[i].x + x_shift;
		y_shifted = inputTrack[i].y + y_shift;

		//Check if in box
		if (abs(x_shifted) < box_edge  && abs(y_shifted) < box_edge)
		{
			//Atomically add to the global counter for the output array length
			outputIndex = atomicAdd(numElementsCompacted,1);
			//printf("current loop i: %d",i);
			//Copy the track inside the box over to the new array
			outputTrack[outputIndex].x = x_shifted;
			outputTrack[outputIndex].y = y_shifted;
			outputTrack[outputIndex].z = inputTrack[i].z;
			outputTrack[outputIndex].edep = inputTrack[i].edep;
		}
	}
}

__global__ void FilterTrackInSphere(double greatestSphereOffset, double sphereRadius, long numSpheresLinear, Track *inputTrack, int *numElements, int *numElementsCompacted, int *trackIdInSphere)
{

	//move all of the variable definitions out of the for loop
	double distFromNearestSphereX, distFromNearestSphereY, distFromNearestSphereZ, dist;

	//Pre-calculate values
	double sphereDiameter = sphereRadius*2;
	double sphereRadiusMag = sphereRadius*sphereRadius; 

	//Determine index and stride
 	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;

	//Loop over all the energy deposition points
	for (long i = index; i < *numElements; i+=stride)
	{		
		//Find the distance from the nearest sphere. You have to shift x_shift by gSO to get in the same coordinate system as the sphere grid
		//An aside: I feel like there is probably a way that you could define the sphere grid that might reduce the complexity of this kernel
		//Another aside: calculating in cubes would reduce complexity as well
		distFromNearestSphereX = llrint((inputTrack[i].x)/sphereDiameter)*sphereDiameter-(inputTrack[i].x);
		distFromNearestSphereY = llrint((inputTrack[i].y)/sphereDiameter)*sphereDiameter-(inputTrack[i].y); 
		distFromNearestSphereZ = llrint((inputTrack[i].z)/sphereDiameter)*sphereDiameter-(inputTrack[i].z); 

		//Determine if inside the nearest sphere
		dist = (distFromNearestSphereX*distFromNearestSphereX)+(distFromNearestSphereY*distFromNearestSphereY)+(distFromNearestSphereZ*distFromNearestSphereZ);

		if (dist <= sphereRadiusMag)
		{
			//Atomically add to the global counter for the output array length
			trackIdInSphere[atomicAdd(numElementsCompacted,1)] = i;
		}
	}
}

__global__ void ScoreTrackInSphere(double greatestSphereOffset, double sphereRadius, long numSpheresLinear, Track *inputTrack, int *numElements, int *trackIdInSphere, long *volumeID, double *edepOutput)
{
	//move all of the variable definitions out of the for loop
	long xIndex, yIndex, zIndex, sphereHitIndex;

	//Pre-calculate some values
	double sphereDiameter = sphereRadius*2;
	double linealDenominator = (2./3.)*sphereDiameter; 

	//Determine index and stride
 	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;

	//Loop over all the energy deposition points
	for (long i = index; i < *numElements; i+=stride)
	{		
		xIndex = llrint((inputTrack[trackIdInSphere[i]].x-greatestSphereOffset)/sphereDiameter);
		yIndex = llrint((inputTrack[trackIdInSphere[i]].y-greatestSphereOffset)/sphereDiameter);
		zIndex = llrint((inputTrack[trackIdInSphere[i]].z-greatestSphereOffset)/sphereDiameter);

		//Determine the Index of the sphere hit
		sphereHitIndex = xIndex + yIndex*(numSpheresLinear) + zIndex*pow(numSpheresLinear,2); //Keep in mind that for the index it starts counting at zero

		//Write to volumeID and edepOutput
		volumeID[i] = sphereHitIndex;
		edepOutput[i] = inputTrack[trackIdInSphere[i]].edep/linealDenominator; //this should be ev/nm which is same a kev/um
	}
}

__global__ void AccumulateHistogramVals(int* temp, int* accumulated,int N)
{
	//Determine index and stride
 	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;

	for (int i = index; i < N; i+=stride)
	{
		accumulated[i] = accumulated[i]+temp[i];
	}
}
__global__ void NonsenseKernel(float *value)
{
	printf("Took in the value:%f", value[0]);
	/*int lastval = 0;
	for (int i = 0; i < 10; i++)
	{
		lastval = atomicAdd(atomicValue,1);
		printf("Last atomic value %d: \n",lastval);
	}*/
}

void testfunction(curandGenerator_t &gen,float *vals)
{
	curandGenerateUniform(gen,vals,2*20000);
	cudaDeviceSynchronize();
}

void mallocfunction(float **vals)
{
	//okay wait, so calling &vals in here, is the pointer to the local object, and not the address of vals
	//Wow! so I pass in a pointer to a pointer. And then feed that to cudaMalloc direclty.
	//Intuitively what this means it memory address pointing to pointer of array of floats --> memory address point to array of floats --> start of floats
	//We had to do this because of the way the malloc works. In c you can't really pass by reference, so when you're mallocing something
	//you need to pass a pointer, to the pointer that you then want to work on
	//I promise it makes sense if you think about it for a moment.

	cudaMalloc(vals,2*sizeof(float)*20000); 
	cudaDeviceSynchronize();
}

void readHostTrack(Track *hostTrack)
{
	for(int i = 0; i < 100000; i++)
	{
		std::cout << hostTrack[i].x << std::endl;
	}
}
__global__ void readDeviceTrack(Track *deviceTrack)
{
	for(int i = 0;i<100000;i++)
	{
		printf("%f \n",deviceTrack[i].x);
	}
}

__global__ void indexTestingKernel()
{
 	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;
	printf("index %d, stride %d \n",index,stride);
}

TH1F score_lineal_GPU(TString filepath, float_t scoring_sphere_spacing, float_t scoring_sphere_diameter, Int_t nthreads, Int_t nSamples = 20000, Long_t random_seed = time(NULL))
{
	//open the file and retrieve the trees
	TFile f = TFile(filepath);
	TTree *trackIndex;
	f.GetObject("Track index",trackIndex);
	long long nTracksInFile = trackIndex->GetEntries();

	//Populate our tuple with the first entry, last entry, and random seed for each thread
	std::vector<std::tuple<Int_t,Int_t,Int_t,TString>> perthread_input_arguments;

	if (nTracksInFile <= nthreads)
	{ 
		long start_entry_val = 0;
		TTreeReader trackIndexReader("Track index", &f);
		TTreeReaderValue<long long> end_entry_val(trackIndexReader, "index");

		for (Int_t i = 0; i < nTracksInFile; i++)
		{
			trackIndexReader.Next();
			perthread_input_arguments.push_back(std::make_tuple(start_entry_val,*end_entry_val-1,i,filepath));
			//Wcout << "thread: " << i << " start val: " << start_entry_val << " end val: " << *end_entry_val-1 << endl;
			start_entry_val = *end_entry_val;
		}
	}
	else
	{
		long start_entry_val = 0;
		TTreeReader trackIndexReader("Track index", &f);
		TTreeReaderValue<long long> end_entry_val(trackIndexReader, "index");

		for (Int_t i = 0; i < nthreads; i++)
		{
			trackIndexReader.Next();
			perthread_input_arguments.push_back(std::make_tuple(start_entry_val,*end_entry_val-1,i,filepath));
			//Wcout << "thread: " << i << " start val: " << start_entry_val << " end val: " << *end_entry_val-1 << endl;
			start_entry_val = *end_entry_val;
		}
		std::cout << "Number of tracks in file greater than requested threads. Case not yet implemented." << std::endl;
	}

	//Pull my geometry information to get passed to each of my threads
	TNamed* voxelSideLength;
	f.GetObject("Voxel side length [mm]",voxelSideLength);
	float scoring_square_half_length = atof(voxelSideLength->GetTitle())*1e6;
	int num_spheres_linear = TMath::Ceil(((scoring_square_half_length*2)/scoring_sphere_spacing)); 
	long long int num_spheres_total = TMath::Power((num_spheres_linear),3);
	double_t top_sphere_offset = -(((float(num_spheres_linear))/2)-0.5)*scoring_sphere_spacing;
	float_t scoringSphereRadius = scoring_sphere_diameter/2;

	//We are done reading the Tree single threaded. Close it.
	f.Close();

	auto workItem = [=](std::tuple<Int_t,Int_t,Int_t,TString> input) //the = sign captures everything in the enclosing function by value. Meaning it makes a process local copy.
	{
		//Calculate size information for memory allocations
		int nVals = std::get<1>(input) - std::get<0>(input) + 1; //+1 because number of values includes first and last value
		size_t trackSize = nVals * sizeof(double); 
		size_t trackStructSize = nVals *sizeof(Track);

		//Define local and GPU track pointers
		Track *hostTrack; Track *deviceTrack; 
		LoadTrack(input, &hostTrack, &deviceTrack); //load track from disk and copy to GPU

		//Allocate memory for the tracks found to be within the box
		Track *inBoxTrack; //In box tracks stores the tracks, with a random shift added, of thos found to be within the box
		cudaMalloc(&inBoxTrack,trackStructSize); 

		//Allocate memory to store the TrackIDs of the points within spheres
		int *inSphereTrackId;
		cudaMalloc(&inSphereTrackId,nVals*sizeof(int));
		
		//Allocate GPU only memory for random numbers
		float *randomVals; 
		GenerateRandomXYShift(input, &randomVals, nSamples, random_seed); //Allocate and fill with random numbers

		//Define histograms
		double *logBins; 
		int* histogramVals; 
		int* histogramValsAccumulated; 
		int nbins = 200; float binLowerMagnitude = -1; float binUpperMagnitude = 2; //Set histogram parameters
		GenerateLogHistogram(&logBins, &histogramVals, &histogramValsAccumulated, nbins, binLowerMagnitude, binUpperMagnitude);
		
		//Allocate GPU only memory for the volume:edep paired list
		long *volumeID; cudaMalloc(&volumeID,sizeof(long)*nVals);
		double *edepInVolume; cudaMalloc(&edepInVolume,trackSize);

		//Allocate memory for the volume:edep paired list after Thrust compacts it
		//TODO: change back to regular memory after debugging
		long *consolidatedVolumeID; cudaMalloc(&consolidatedVolumeID,sizeof(long)*nVals);
		double *consolidatedEdepInVolume; cudaMalloc(&consolidatedEdepInVolume,trackSize);
		
		//Todo: Create the output vectors for after thrust has summed my raw data into a consolidated list of volumeID:edep, rather than redifining each time

		int *NumInBox; int *NumInSpheres;
		cudaMallocManaged(&NumInBox,sizeof(int)); cudaMallocManaged(&NumInSpheres,sizeof(int));

		for (int j = 0; j < nSamples; j++)
		{
			//Allocate histogram memory
			void* histogramTempStorage = NULL;
			size_t tempStorageSize = 0;

			//New track. Zero values
			cudaMemset(NumInBox,0,sizeof(int));
			cudaMemset(NumInSpheres,0,sizeof(int));

			//indexTestingKernel<<<32,32>>>();

			FilterInScoringBox<<<60,256>>>(top_sphere_offset,scoringSphereRadius,num_spheres_linear,randomVals,deviceTrack,inBoxTrack,nVals,NumInBox,j);	
			FilterTrackInSphere<<<60,256>>>(top_sphere_offset,scoringSphereRadius,num_spheres_linear,inBoxTrack,NumInBox,NumInSpheres,inSphereTrackId);
			ScoreTrackInSphere<<<60,256>>>(top_sphere_offset,scoringSphereRadius,num_spheres_linear,inBoxTrack,NumInSpheres,inSphereTrackId,volumeID,edepInVolume);
			
			//I think, because the sort_by_key operation takes *NumInSpheres as an argument
			//If the kernel call is given, before NumInSpheres has finished updating, then it gets an incorrect value
			cudaDeviceSynchronize();

			//Use Thrust, to sort my energy depositions in the order of the volumes they occured in 
			thrust::sort_by_key(thrust::device,volumeID,volumeID+*NumInSpheres,edepInVolume);
			thrust::pair<long*,double*> endOfReducedList;
			endOfReducedList = thrust::reduce_by_key(thrust::device,volumeID,volumeID+*NumInSpheres,edepInVolume,consolidatedVolumeID,consolidatedEdepInVolume); //Then reduce the energy depositions. Default reduction function is plus(), which is exactly what I want. i.e. summing the depositions
			
			//First call to the histogram allocates the temp storage and size
			cub::DeviceHistogram::HistogramRange(histogramTempStorage,tempStorageSize,consolidatedEdepInVolume,histogramVals,nbins+1,logBins,endOfReducedList.second-consolidatedEdepInVolume);

			//Allocate the temporary storage
			cudaMalloc(&histogramTempStorage,tempStorageSize); 

			//Second call populates the histogram
			cub::DeviceHistogram::HistogramRange(histogramTempStorage,tempStorageSize,consolidatedEdepInVolume,histogramVals,nbins+1,logBins,endOfReducedList.second-consolidatedEdepInVolume);

			//Accumulate the histogram values
			AccumulateHistogramVals<<<4,32>>>(histogramVals,histogramValsAccumulated,nbins);
			
			//std::cout << *NumInBox << " " << *NumInSpheres << std::endl;
		}

		int number_of_values_in_histogram = 0;
		cudaDeviceSynchronize();
		//Read out histogram
		for (int i = 0; i < nbins; i++)
		{
			number_of_values_in_histogram += histogramValsAccumulated[i];
			std::cout << "Bin: " << logBins[i] << " Counts: " << histogramValsAccumulated[i] << std::endl;
		}

		std::cout << number_of_values_in_histogram << std::endl;
		//TODO: close my file at some point

	  	//Initialize the histogram
		TH1F lineal_histogram = TH1F("Lineal energy histogram", "y*f(y)", 200, -2,1);
		return lineal_histogram;

	};


	// Create the pool of workers
	ROOT::TProcessExecutor workers(nthreads);
	//Process the jobs and get a vector of the output
	std::vector<TH1F> process_output = workers.Map(workItem, perthread_input_arguments);

	TH1F lineal_histogram = TH1F("Lineal energy histogram", "y*f(y)", 200, -2,1);

	return lineal_histogram;

}

//TODO: Change this to work with a C-style struct later, so x,y,z,edep are all one entry
/*__global__ void SuperimposeTrack(double greatestSphereOffset, double sphereDiameter, long numSpheresLinear, float* randomVals, double* x, double* y, double* z, double* edep,long *volumeID, double *edepOutput, long numElements,int oversampleIterationNumber)
{
	//Our entire geometry should be able to be described by only the greatest offset, the sphere diameter and number of spheres in a line. That's useful
	double sphereRadius = sphereDiameter/2;
	double linealDenominator = (2./3.)*sphereDiameter; //calculate this here as an efficiency gain

	//Determine index and stride
 	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;

	//Convert random shifts  in to appropriate range
	double x_shift = ((randomVals[(oversampleIterationNumber*2)]*greatestSphereOffset*2)-greatestSphereOffset);
	double y_shift = ((randomVals[(oversampleIterationNumber*2+1)]*greatestSphereOffset*2)-greatestSphereOffset);
	//printf("x_shift: %f \n",randomVals[(oversampleIterationNumber*2)]);
	//printf("y_shift: %f \n",randomVals[(oversampleIterationNumber*2+1)]);
	//printf("Greatest sphere offset: %f \n", greatestSphereOffset);

	//Loop over all the energy deposition points
	for (long i = index; i < numElements; i+=stride)
	{
		//Write a zero to edepOutput and volumeID. Doing this here avoids warp divergence later.
		edepOutput[i] = 0; volumeID[i] = 0;

		//Apply random shift. My numbers that come in are floats from 0.0 to 1.0. Have to shift them to the desired range
		double x_shifted = x[i] + x_shift;
		double y_shifted = y[i] + y_shift;

		//Check if inside box
		if (abs(x_shifted) < abs(greatestSphereOffset)+(sphereRadius) && abs(y_shifted) < abs(greatestSphereOffset)+(sphereRadius) && abs(z[i]) < abs(greatestSphereOffset)+(sphereRadius))
		{
			//Convert position to index in the grid of spheres
			//printf("x_shifted: %f \n",x_shifted);
			long xIndex = llround((x_shifted-greatestSphereOffset)/sphereDiameter);
			long yIndex = llround((y_shifted-greatestSphereOffset)/sphereDiameter);
			long zIndex = llround((z[i]-greatestSphereOffset)/sphereDiameter);
			
			//Determine the location of the nearest sphere in the grid (with 0,0,0 being the top left sphere, different coordinate system than the ptcls are in)
			double nearestSphereX = xIndex*sphereDiameter;
			double nearestSphereY = yIndex*sphereDiameter;
			double nearestSphereZ = zIndex*sphereDiameter;

			//Find the distance from the nearest sphere. You have to shift x_shift by gSO to get in the same coordinate system as the sphere grid
			//An aside: I feel like there is probably a way that you could define the sphere grid that might reduce the complexity of this kernel
			//Another aside: calculating in cubes would reduce complexity as well
			double distFromNearestSphereX = nearestSphereX-(x_shifted-greatestSphereOffset);
			double distFromNearestSphereY = nearestSphereY-(y_shifted-greatestSphereOffset); 
			double distFromNearestSphereZ = nearestSphereZ-(z[i]-greatestSphereOffset); 

			//Determine if inside the nearest sphere
			double dist = pow(distFromNearestSphereX,2)+pow(distFromNearestSphereY,2)+pow(distFromNearestSphereZ,2);
			dist = sqrt(dist);

			if (dist <= sphereRadius)
			{
				//Determine the Index of the sphere hit
				long sphereHitIndex = xIndex + yIndex*(numSpheresLinear) + zIndex*pow(numSpheresLinear,2); //Keep in mind that for the index it starts counting at zero

				//Write to volumeID and edepOutput
				volumeID[i] = sphereHitIndex;
				edepOutput[i] = edep[i]/linealDenominator; //this should be ev/nm which is same a kev/um
				//printf("volumeID: %ld %ld %ld edep: %f. \n",xIndex,yIndex,zIndex,edepOutput[i]);
			}
		}
	}
}*/

/*__global__ void ScoreTrackInSphere(double greatestSphereOffset, double sphereRadius, long numSpheresLinear, Track *inputTrack, int *numElements, long *volumeID, double *edepOutput, int *numElementsCompacted)
{

	double sphereDiameter = sphereRadius*2;
	double linealDenominator = (2./3.)*sphereDiameter; 

	//Determine index and stride
 	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;

	//move all of the variable definitions out of the for loop
	long xIndex, yIndex, zIndex, sphereHitIndex;
	double nearestSphereX, nearestSphereY, nearestSphereZ, distFromNearestSphereX, distFromNearestSphereY, distFromNearestSphereZ, dist;
	double xRelativeToEdge, yRelativeToEdge, zRelativeToEdge;
	int outputIndex;

	//Loop over all the energy deposition points
	for (long i = index; i < *numElements; i+=stride)
	{
		xRelativeToEdge = inputTrack[i].x-greatestSphereOffset;
		yRelativeToEdge = inputTrack[i].y-greatestSphereOffset;
		zRelativeToEdge = inputTrack[i].z-greatestSphereOffset;

		//Convert position to index in the grid of spheres
		//printf("x_shifted: %f \n",x_shifted);
		xIndex = llround((xRelativeToEdge)/sphereDiameter);
		yIndex = llround((yRelativeToEdge)/sphereDiameter);
		zIndex = llround((zRelativeToEdge)/sphereDiameter);
		
		//Determine the location of the nearest sphere in the grid (with 0,0,0 being the top left sphere, different coordinate system than the ptcls are in)
		nearestSphereX = xIndex*sphereDiameter;
		nearestSphereY = yIndex*sphereDiameter;
		nearestSphereZ = zIndex*sphereDiameter;

		//Find the distance from the nearest sphere. You have to shift x_shift by gSO to get in the same coordinate system as the sphere grid
		//An aside: I feel like there is probably a way that you could define the sphere grid that might reduce the complexity of this kernel
		//Another aside: calculating in cubes would reduce complexity as well
		distFromNearestSphereX = nearestSphereX-(xRelativeToEdge);
		distFromNearestSphereY = nearestSphereY-(yRelativeToEdge); 
		distFromNearestSphereZ = nearestSphereZ-(zRelativeToEdge); 

		//Determine if inside the nearest sphere
		dist = pow(distFromNearestSphereX,2)+pow(distFromNearestSphereY,2)+pow(distFromNearestSphereZ,2);
		dist = sqrt(dist);

		if (dist <= sphereRadius)
		{
			//Determine the Index of the sphere hit
			sphereHitIndex = xIndex + yIndex*(numSpheresLinear) + zIndex*pow(numSpheresLinear,2); //Keep in mind that for the index it starts counting at zero
			
			//Atomically add to the global counter for the output array length
			outputIndex = atomicAdd(numElementsCompacted,1);

			//Write to volumeID and edepOutput
			volumeID[outputIndex] = sphereHitIndex;
			edepOutput[outputIndex] = inputTrack[i].edep/linealDenominator; //this should be ev/nm which is same a kev/um
			//printf("volumeID: %ld %ld %ld edep: %f. \n",xIndex,yIndex,zIndex,edepOutput[i]);
		}
	}
}*/