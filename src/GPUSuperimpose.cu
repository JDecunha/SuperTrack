//SuperTrack
#include "testCUDA.cuh"
#include "utils.hh"
#include "SuperTrackTypes.cuh"
//ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TEntryList.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TMath.h"
#include "ROOT/TProcessExecutor.hxx"
//STD
#include <vector>
#include <iterator>
#include <tuple>
#include <filesystem>
//CUDA Libraries
#include <cuda.h>
#include <curand.h>
//CUB (Cuda UnBound)
#include <cub/cub.cuh>


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

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

__global__ void FilterInScoringBox(SphericalGeometry geometry, float* randomVals, Track *inputTrack, Track *outputTrack, int numElements, int *numElementsCompacted, int oversampleIterationNumber)
{
	//This function
	//1.) Applies the random shift to the x,y coordinates
	//2.) Checks which edep events are in the box
	//3.) Performs stream compaction on those events which are in the box
	//We are using stream compaction to avoid a monolithic kernel with large blocks within if-statements which reduces warp efficiency

	//Put definitions outside of for-loop to prevent repeat constructor calls
	double x_shifted; double y_shifted; int outputIndex;

	//Determine index and strid
  int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;

	//Counters for shared memory atomics
	int localPosition;
	__shared__ int localIndexCounter;
	
	//Convert random shifts in to appropriate range
	double x_shift = ((randomVals[(oversampleIterationNumber*2)]*geometry.greatestSphereOffset*2)-geometry.greatestSphereOffset);
	double y_shift = ((randomVals[(oversampleIterationNumber*2+1)]*geometry.greatestSphereOffset*2)-geometry.greatestSphereOffset);

	//The value we compare to, to check if it's in the box
	double box_edge = abs(geometry.greatestSphereOffset)+(geometry.sphereRadius);

	//Loop over all the energy deposition points
	for (int i = index; i < numElements; i+=stride)
	{
		//Apply random shift
		x_shifted = inputTrack[i].x + x_shift;
		y_shifted = inputTrack[i].y + y_shift;

		//Set local position to negative value, only takes on positive value if predicate is true
		localPosition = -1;

		//Zero the local counter
		if (threadIdx.x == 0) 
		{
			localIndexCounter = 0;
		}
		__syncthreads();

		//Check if in box, if true assign the local index position
		if (abs(x_shifted) < box_edge  && abs(y_shifted) < box_edge)
		{
			localPosition = atomicAdd(&localIndexCounter,1);
		}
		__syncthreads();

		//Add the local counter to the global counter
		if (threadIdx.x == 0)
		{
			localIndexCounter = atomicAdd(numElementsCompacted,localIndexCounter);
		}
		__syncthreads();

		//If predicate is true, then write the track to position localCounter+localPosition (localCounter now stores the globalCounter value because of the atomic add)
		if(localPosition != -1)
		{
			//Atomically add to the global counter for the output array length
			outputIndex = localPosition+localIndexCounter;
			//printf("current loop i: %d",i);
			//Copy the track inside the box over to the new array
			outputTrack[outputIndex].x = x_shifted;
			outputTrack[outputIndex].y = y_shifted;
			outputTrack[outputIndex].z = inputTrack[i].z;
			outputTrack[outputIndex].edep = inputTrack[i].edep;
		}
		__syncthreads();
	}
}

__global__ void FilterTrackInSphere(SphericalGeometry geometry, Track *inputTrack, int *numElements, int *numElementsCompacted, int *trackIdInSphere)
{

	//printf("%f %f \n",geometry.sphereDiameter,geometry.scoringRegionHalfLength);

	//move all of the variable definitions out of the for loop
	double distFromNearestSphereX, distFromNearestSphereY, distFromNearestSphereZ, dist;

	//Determine index and stride
 	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;

	//Counters for shared memory atomics
	int localPosition;
	__shared__ int localIndexCounter;

	//Pre-calculate values
	double sphereDiameter = geometry.sphereDiameter; 
	double sphereRadiusMag = geometry.sphereRadius*geometry.sphereRadius; 

	//Loop over all the energy deposition points
	for (long i = index; i < *numElements; i+=stride)
	{		
		//Find distance to the nearest sphere.
		//For performance reasons, we work in an arbitrary coordinate system here, rather than the "global" coordinate system
		//The "global" coordinate systrem is relative to the greatest sphere offset
		//In the later scoring kernel we work in the global coordinate system, and that's why we subtract the greatest sphere offset there
		distFromNearestSphereX = llrint((inputTrack[i].x)/sphereDiameter)*geometry.sphereDiameter-(inputTrack[i].x);
		distFromNearestSphereY = llrint((inputTrack[i].y)/sphereDiameter)*geometry.sphereDiameter-(inputTrack[i].y); 
		distFromNearestSphereZ = llrint((inputTrack[i].z)/sphereDiameter)*geometry.sphereDiameter-(inputTrack[i].z); 

		//Determine if inside the nearest sphere
		dist = (distFromNearestSphereX*distFromNearestSphereX)+(distFromNearestSphereY*distFromNearestSphereY)+(distFromNearestSphereZ*distFromNearestSphereZ);

		//Set local position to negative value, only takes on positive value if predicate is true
		localPosition = -1;

		//Zero the local counter
		if (threadIdx.x == 0) 
		{
			localIndexCounter = 0;
		}
		__syncthreads();

		//Check if in sphere, then assign local index position
		if (dist <= sphereRadiusMag)
		{
			localPosition = atomicAdd(&localIndexCounter,1);
			
		}
		__syncthreads();

		//Add the local counter to the global counter
		if (threadIdx.x == 0)
		{
			localIndexCounter = atomicAdd(numElementsCompacted,localIndexCounter);
		}
		__syncthreads();

		//If predicate is true, then write the track to position localCounter+localPosition (localCounter now stores the globalCounter value because of the atomic add)
		if (localPosition != -1)
		{
			//Atomically add to the global counter for the output array length
			trackIdInSphere[localPosition+localIndexCounter] = i;
		}
		__syncthreads();

	}
}

__global__ void ScoreTrackInSphere(SphericalGeometry geometry, Track *inputTrack, int *numElements, int *trackIdInSphere, VolumeEdepPair outputPair)
{
	//move all of the variable definitions out of the for loop
	long xIndex, yIndex, zIndex, sphereHitIndex;

	//Pre-calculate some values
	double sphereDiameter = geometry.sphereDiameter;
	double linealDenominator = (2./3.)*sphereDiameter; 

	//Determine index and stride
 	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;

	//Loop over all the energy deposition points
	for (uint64_t i = index; i < *numElements; i+=stride)
	{		
		xIndex = llrint((inputTrack[trackIdInSphere[i]].x-geometry.greatestSphereOffset)/sphereDiameter);
		yIndex = llrint((inputTrack[trackIdInSphere[i]].y-geometry.greatestSphereOffset)/sphereDiameter);
		zIndex = llrint((inputTrack[trackIdInSphere[i]].z-geometry.greatestSphereOffset)/sphereDiameter);

		//Determine the Index of the sphere hit
		sphereHitIndex = xIndex + yIndex*(geometry.numSpheresLinear) + zIndex*pow(geometry.numSpheresLinear,2); //Keep in mind that for the index it starts counting at zero

		//Write to volumeID and edepOutput
		outputPair.volume[i] = sphereHitIndex;
		outputPair.edep[i] = inputTrack[trackIdInSphere[i]].edep/linealDenominator; //this should be ev/nm which is same a kev/um
	}
}

struct CUBAddOperator
{
    template <typename T>
    CUB_RUNTIME_FUNCTION __forceinline__
    T operator()(const T &a, const T &b) const {
        return a+b;
    }
};

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
__global__ void readDeviceEdepList(double *edeps)
{
	for(int i = 0;i<200;i++)
	{
		printf("%f \n",edeps[i]);
	}
}

__global__ void indexTestingKernel()
{
 	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;
	printf("index %d, stride %d \n",index,stride);
}

/*__global__ void extractCudaInformation()
{
	printf("Cuda_ARCH: %s",__CUDA_ARCH__);
	printf("CUDA_ACC: %s",__CUDACC__);
	printf("CUDA_RDC: %s",__CUDACC_RDC__);
}*/

//TODO:Change this back to Malloc after testing
VolumeEdepPair AllocateGPUVolumeEdepPair(uint64_t numElements)
{
	VolumeEdepPair toAllocate;
	cudaMalloc(&(toAllocate.volume),numElements*sizeof(uint64_t));
	cudaMalloc(&(toAllocate.edep),numElements*sizeof(double));
	cudaMallocManaged(&(toAllocate.numElements),sizeof(int));

	return toAllocate;
}

//Todo: create function for freeing GPUVolumeEdepPair

CubStorageBuffer AllocateCubSortBuffer(VolumeEdepPair edepPairList, uint64_t nVals)
{
	//Create the buffer with default constructor
	CubStorageBuffer returnBuffer = CubStorageBuffer();

	//Call the CUB function to determine memory constraints, then malloc
	cub::DeviceRadixSort::SortPairs(returnBuffer.storage,returnBuffer.size,edepPairList.volume,edepPairList.volume,edepPairList.edep,edepPairList.edep,nVals);
	cudaMalloc(&returnBuffer.storage,returnBuffer.size); 

	return returnBuffer;
}

CubStorageBuffer AllocateCubReduceBuffer(VolumeEdepPair edepPairList, uint64_t nVals)
{
	//Create the buffer with default constructor
	CubStorageBuffer returnBuffer = CubStorageBuffer();
	//Generic Reduction Operator
	CUBAddOperator reductionOperator;

	//Call the CUB function to determine memory constraints, then malloc
	cub::DeviceReduce::ReduceByKey(returnBuffer.storage,returnBuffer.size, edepPairList.volume, edepPairList.volume, edepPairList.edep, edepPairList.edep, edepPairList.numElements, reductionOperator, nVals);
	cudaMalloc(&returnBuffer.storage,returnBuffer.size); 

	return returnBuffer;
}

CubStorageBuffer AllocateCubHistogramBuffer(VolumeEdepPair edepPairList, uint64_t nVals, int* histogramVals, double* logBins, int nbins)
{
	//Create the buffer with default constructor
	CubStorageBuffer returnBuffer = CubStorageBuffer();

	//Call the CUB function to determine memory constraints, then malloc
	cub::DeviceHistogram::HistogramRange(returnBuffer.storage,returnBuffer.size, edepPairList.edep,histogramVals,nbins+1,logBins,nVals);
	cudaMalloc(&returnBuffer.storage,returnBuffer.size); 

	return returnBuffer;
}

__global__ void ReadVolumeEdepPair(VolumeEdepPair* pair)
{
	for (int i = 0; i < 10000; i++)
	{
		printf("volume: %d",pair->volume[i]);
	}	
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

	SphericalGeometry sphericalGeometry = SphericalGeometry(scoring_square_half_length,scoring_sphere_diameter);

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
		Track *randomlyShiftedTrack; 
		cudaMalloc(&randomlyShiftedTrack,trackStructSize); 

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
		VolumeEdepPair edepsInTarget = AllocateGPUVolumeEdepPair(nVals);
		VolumeEdepPair sortedEdeps = AllocateGPUVolumeEdepPair(nVals); 
		VolumeEdepPair reducedEdeps = AllocateGPUVolumeEdepPair(nVals); 

		int *NumInBox; 
		cudaMallocManaged(&NumInBox,sizeof(int)); 

		CUBAddOperator reductionOperator;

		//Allocate memory for the temporary storage the CUB operations needs
		CubStorageBuffer sortBuffer = AllocateCubSortBuffer(edepsInTarget,nVals);
		CubStorageBuffer reduceBuffer = AllocateCubReduceBuffer(edepsInTarget,nVals);
		CubStorageBuffer histogramBuffer = AllocateCubHistogramBuffer(edepsInTarget,nVals,histogramVals,logBins,nbins);


		for (int j = 0; j < nSamples; j++)
		{

			//New track. Zero values
			cudaMemset(NumInBox,0,sizeof(int));
			cudaMemset(edepsInTarget.numElements,0,sizeof(int));

			FilterInScoringBox<<<60,256>>>(sphericalGeometry,randomVals,deviceTrack,randomlyShiftedTrack,nVals,NumInBox,j);	
			FilterTrackInSphere<<<60,256>>>(sphericalGeometry,randomlyShiftedTrack,NumInBox,edepsInTarget.numElements,inSphereTrackId); 
			ScoreTrackInSphere<<<60,256>>>(sphericalGeometry,randomlyShiftedTrack,edepsInTarget.numElements,inSphereTrackId,edepsInTarget); 
		
			cudaDeviceSynchronize();

			//I think, because the sort_by_key operation takes *NumInSpheres as an argument
			//If the kernel call is given, before NumInSpheres has finished updating, then it gets an incorrect value

			//Sort the edep volume pairs
			cub::DeviceRadixSort::SortPairs(sortBuffer.storage,sortBuffer.size,edepsInTarget.volume,sortedEdeps.volume,edepsInTarget.edep,sortedEdeps.edep,*(edepsInTarget.numElements));

			// reduce the energy depositions
			cub::DeviceReduce::ReduceByKey(reduceBuffer.storage,reduceBuffer.size, sortedEdeps.volume, reducedEdeps.volume, sortedEdeps.edep, reducedEdeps.edep, reducedEdeps.numElements, reductionOperator, *(edepsInTarget.numElements));

			//Create the histogram
			cub::DeviceHistogram::HistogramRange(histogramBuffer.storage,histogramBuffer.size,reducedEdeps.edep,histogramVals,nbins+1,logBins,*reducedEdeps.numElements);

			//Accumulate the histogram values
			AccumulateHistogramVals<<<4,32>>>(histogramVals,histogramValsAccumulated,nbins);
	
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
		//TODO: close my file at some point */
		//cudaDeviceSynchronize();
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