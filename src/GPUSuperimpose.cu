//SuperTrack
#include "utils.hh"
#include "SuperTrackTypes.cuh"
#include "Track.cuh"
#include "CubStorageBuffer.cuh"
#include "Histogram.cuh"
#include "VolumeEdepPair.cuh"
#include "ThreadAllocation.hh"
//#include "SphericalGeometryWithBounding.cuh"
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

__global__ void FilterInScoringBox(SphericalGeometry geometry, float* randomVals, Track inputTrack, Track outputTrack, int numElements, int *numElementsCompacted, int oversampleIterationNumber)
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
		x_shifted = inputTrack.x[i] + x_shift;
		y_shifted = inputTrack.y[i] + y_shift;

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
			outputTrack.x[outputIndex] = x_shifted;
			outputTrack.y[outputIndex] = y_shifted;
			outputTrack.z[outputIndex] = inputTrack.z[i];
			outputTrack.edep[outputIndex] = inputTrack.edep[i];
		}
		__syncthreads();
	}
}

__global__ void FilterTrackInSphere(SphericalGeometry geometry, Track inputTrack, int *numElements, int *numElementsCompacted, int *trackIdInSphere)
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
		distFromNearestSphereX = llrint((inputTrack.x[i])/sphereDiameter)*geometry.sphereDiameter-(inputTrack.x[i]);
		distFromNearestSphereY = llrint((inputTrack.y[i])/sphereDiameter)*geometry.sphereDiameter-(inputTrack.y[i]); 
		distFromNearestSphereZ = llrint((inputTrack.z[i])/sphereDiameter)*geometry.sphereDiameter-(inputTrack.z[i]); 

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

__global__ void ScoreTrackInSphere(SphericalGeometry geometry, Track inputTrack, int *numElements, int *trackIdInSphere, VolumeEdepPair outputPair)
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
		xIndex = llrint((inputTrack.x[trackIdInSphere[i]]-geometry.greatestSphereOffset)/sphereDiameter);
		yIndex = llrint((inputTrack.y[trackIdInSphere[i]]-geometry.greatestSphereOffset)/sphereDiameter);
		zIndex = llrint((inputTrack.z[trackIdInSphere[i]]-geometry.greatestSphereOffset)/sphereDiameter);

		//Determine the Index of the sphere hit
		sphereHitIndex = xIndex + yIndex*geometry.numSpheresLinear+ zIndex*geometry.numSpheresLinear*geometry.numSpheresLinear; //Keep in mind that for the index it starts counting at zero

		//Write to volumeID and edepOutput
		outputPair.volume[i] = sphereHitIndex;
		outputPair.edep[i] = inputTrack.edep[trackIdInSphere[i]]/linealDenominator; //this should be ev/nm which is same a kev/um
	}
}

__global__ void ZeroInt(int* toZero)
{
	*toZero = 0;
}

TH1F score_lineal_GPU_New(std::vector<ThreadAllocation> threadAllocations, float_t scoring_sphere_spacing, float_t scoring_sphere_diameter)
{
	//open the file and retrieve the trees
	TFile f = TFile((threadAllocations[0].GetTasks())[0].GetFilename());

	//Pull my geometry information to get passed to each of my threads
	TNamed* voxelSideLength;
	f.GetObject("Voxel side length [mm]",voxelSideLength);
	float scoring_square_half_length = atof(voxelSideLength->GetTitle())*1e6;
	SphericalGeometry sphericalGeometry = SphericalGeometry(scoring_square_half_length,scoring_sphere_diameter);

	//We are done reading the Tree single threaded. Close it.
	f.Close();

	auto workItem = [=](ThreadAllocation threadInput) //the = sign captures everything in the enclosing function by value. Meaning it makes a process local copy.
	{
		std::vector<ThreadTask> tasks = threadInput.GetTasks();

		//New design idea: Loop over an iterator i
		//Then each of my functions can take (ThreadAllocation allocation, int taskNumber)
		//That's how I'll pass all the information I need to these classes

		for (ThreadTask task : tasks) //loop over every track requested
		{
			int nVals = (task.GetExitPoint() - task.GetEntryPoint()) + 1; //+1 because number of values includes first and last value

			std::tuple<Int_t,Int_t,Int_t,TString> input = std::make_tuple(task.GetEntryPoint(),task.GetExitPoint(),threadInput.GetRandomSeed(),task.GetFilename());

			//Define Track
			Track deviceTrack;
			deviceTrack.AllocateAndLoadTrack(input); //load track from disk and copy to GPU

			//Allocate memory for the tracks found to be within the box
			Track randomlyShiftedTrack;
			randomlyShiftedTrack.AllocateEmptyTrack(nVals);

			//Allocate memory to store the TrackIDs of the points within spheres
			int *inSphereTrackId;
			cudaMalloc(&inSphereTrackId,nVals*sizeof(int));
			
			//Allocate GPU only memory for random numbers
			float *randomVals; 
			GenerateRandomXYShift(input, &randomVals, threadInput.GetNOversamples(), threadInput.GetRandomSeed()); //Allocate and fill with random numbers

			//Allocate GPU only memory for the volume:edep paired list
			VolumeEdepPair edepsInTarget;
			edepsInTarget.Allocate(nVals);

			//Make histogram
			Histogram histogram = Histogram(200,-1,2,"log");
			histogram.Allocate(edepsInTarget);
			
			int *NumInBox; 
			cudaMallocManaged(&NumInBox,sizeof(int)); 

			//Configure cuda kernel launches
			/*int blockSize;
			int minGridSize;
			int gridSize;
			cudaOccupancyMaxPotentialBlockSize(&minGridSize,&blockSize,ScoreTrackInSphere,0,0);
			gridSize = (nVals + blockSize - 1)/blockSize;
			std::cout << gridSize << " " << blockSize << std::endl;*/

			for (int j = 0; j < threadInput.GetNOversamples(); j++)
			{
				//New track. Zero values
				ZeroInt<<<1,1>>>(NumInBox);
				ZeroInt<<<1,1>>>(edepsInTarget.numElements);

				//Filter and score the tracks
				FilterInScoringBox<<<256,256>>>(sphericalGeometry,randomVals,deviceTrack,randomlyShiftedTrack,nVals,NumInBox,j);	
				FilterTrackInSphere<<<256,256>>>(sphericalGeometry,randomlyShiftedTrack,NumInBox,edepsInTarget.numElements,inSphereTrackId); 
				ScoreTrackInSphere<<<256,256>>>(sphericalGeometry,randomlyShiftedTrack,edepsInTarget.numElements,inSphereTrackId,edepsInTarget); 

				//Sort the edeps by volumeID, reduce (accumulate), and then place into histograms
				histogram.SortReduceAndAddToHistogram(edepsInTarget);
			}

			int number_of_values_in_histogram = 0;
			cudaDeviceSynchronize();
			//Read out histogram
			for (int i = 0; i < histogram._nbins; i++)
			{
				number_of_values_in_histogram += histogram._histogramValsAccumulated[i];
				std::cout << "Bin: " << histogram._binEdges[i] << " Counts: " << histogram._histogramValsAccumulated[i] << std::endl;
			}

			std::cout << number_of_values_in_histogram << std::endl;

			//Free directly allocated memory
			cudaFree(inSphereTrackId);
			cudaFree(randomVals);
			cudaFree(NumInBox);

			//Free my classes
			edepsInTarget.Free();
			deviceTrack.Free();
			randomlyShiftedTrack.Free();
			histogram.Free();

			//Initialize the histogram
			TH1F lineal_histogram = TH1F("Lineal energy histogram", "y*f(y)", 200, -2,1);
			return lineal_histogram;
		}


	};

	// Create the pool of workers
	ROOT::TProcessExecutor workers(threadAllocations.size());
	//Process the jobs and get a vector of the output
	std::vector<TH1F> process_output = workers.Map(workItem, threadAllocations);

	TH1F lineal_histogram = TH1F("Lineal energy histogram", "y*f(y)", 200, -2,1);

	return lineal_histogram;

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
	SphericalGeometry sphericalGeometry = SphericalGeometry(scoring_square_half_length,scoring_sphere_diameter);

	//We are done reading the Tree single threaded. Close it.
	f.Close();

	auto workItem = [=](std::tuple<Int_t,Int_t,Int_t,TString> input) //the = sign captures everything in the enclosing function by value. Meaning it makes a process local copy.
	{
		int nVals = std::get<1>(input) - std::get<0>(input) + 1; //+1 because number of values includes first and last value

		//Define Track
		Track deviceTrack;
		deviceTrack.AllocateAndLoadTrack(input); //load track from disk and copy to GPU

		//Allocate memory for the tracks found to be within the box
		Track randomlyShiftedTrack;
		randomlyShiftedTrack.AllocateEmptyTrack(nVals);

		//Allocate memory to store the TrackIDs of the points within spheres
		int *inSphereTrackId;
		cudaMalloc(&inSphereTrackId,nVals*sizeof(int));
		
		//Allocate GPU only memory for random numbers
		float *randomVals; 
		GenerateRandomXYShift(input, &randomVals, nSamples, random_seed); //Allocate and fill with random numbers

		//Allocate GPU only memory for the volume:edep paired list
		VolumeEdepPair edepsInTarget;
		edepsInTarget.Allocate(nVals);

		//Make histogram
		Histogram histogram = Histogram(200,-1,2,"log");
		histogram.Allocate(edepsInTarget);
		
		int *NumInBox; 
		cudaMallocManaged(&NumInBox,sizeof(int)); 

		//Configure cuda kernel launches
		/*int blockSize;
		int minGridSize;
		int gridSize;
		cudaOccupancyMaxPotentialBlockSize(&minGridSize,&blockSize,ScoreTrackInSphere,0,0);
		gridSize = (nVals + blockSize - 1)/blockSize;
		std::cout << gridSize << " " << blockSize << std::endl;*/

		for (int j = 0; j < nSamples; j++)
		{
			//New track. Zero values
			ZeroInt<<<1,1>>>(NumInBox);
			ZeroInt<<<1,1>>>(edepsInTarget.numElements);

			//Filter and score the tracks
			FilterInScoringBox<<<256,256>>>(sphericalGeometry,randomVals,deviceTrack,randomlyShiftedTrack,nVals,NumInBox,j);	
			FilterTrackInSphere<<<256,256>>>(sphericalGeometry,randomlyShiftedTrack,NumInBox,edepsInTarget.numElements,inSphereTrackId); 
			ScoreTrackInSphere<<<256,256>>>(sphericalGeometry,randomlyShiftedTrack,edepsInTarget.numElements,inSphereTrackId,edepsInTarget); 

			//Sort the edeps by volumeID, reduce (accumulate), and then place into histograms
			histogram.SortReduceAndAddToHistogram(edepsInTarget);
		}

		int number_of_values_in_histogram = 0;
		cudaDeviceSynchronize();
		//Read out histogram
		for (int i = 0; i < histogram._nbins; i++)
		{
			number_of_values_in_histogram += histogram._histogramValsAccumulated[i];
			std::cout << "Bin: " << histogram._binEdges[i] << " Counts: " << histogram._histogramValsAccumulated[i] << std::endl;
		}

		std::cout << number_of_values_in_histogram << std::endl;
		//TODO: close my file at some point */
	  


		//Free directly allocated memory
		cudaFree(inSphereTrackId);
		cudaFree(randomVals);
		cudaFree(NumInBox);

		//Free my classes
		edepsInTarget.Free();
		deviceTrack.Free();
		randomlyShiftedTrack.Free();
		histogram.Free();

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