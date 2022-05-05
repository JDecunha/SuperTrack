#include "VoxelConstrainedSphereMethod.hh"
#include "VolumeEdepPair.cuh"
#include "Track.cuh"
#include <cuda.h>
#include <curand.h>
#include <iostream>

//Constructor
VoxelConstrainedSphereMethod::VoxelConstrainedSphereMethod(const INIReader& macroReader) : SimulationMethod(macroReader), _sphericalGeometry(macroReader)
{
	ParseInput();
}

//ParseInput takes the INIReader to initialize the class
void VoxelConstrainedSphereMethod::ParseInput() //Many of the inputs are currently handled by the SphericalGeometry helper struct, look there as well
{
	//Pull the number of SMs 
	int deviceId;
	cudaGetDevice(&deviceId);
	cudaDeviceProp props;
	cudaGetDeviceProperties(&props, deviceId);

	//Set default number of Blocks and Threads
	int defaultNumBlocks = props.multiProcessorCount*10;
	int defaultNumThreads = 256;

	_suggestedCudaBlocks = _macroReader.GetReal("VoxelConstrainedSphere","SuggestedCudaBlocks",defaultNumBlocks);
	_suggestedCudaThreads = _macroReader.GetReal("VoxelConstrainedSphere","SuggestedCudaThreads",defaultNumThreads);
	_randomShifts = _macroReader.GetBoolean("VoxelConstrainedSphere","ShiftTrack",true);
}

//Static method that the SimulationMethodFactory uses to build this simulation method
SimulationMethod* VoxelConstrainedSphereMethod::Construct(const INIReader& macroReader)
{
	return new VoxelConstrainedSphereMethod(macroReader);
}

//Called to allocate memory at the start of processing a track
void VoxelConstrainedSphereMethod::AllocateTrackProcess(Track track, ThreadTask task) 
{ 
	_oversampleIterationNumber = 0;
	_nSteps = task.GetExitPoint() - task.GetEntryPoint();

	//Allocate GPU only memory and fill with random numbers, or fill with 0.5 if shifts not requested
	GenerateRandomXYShift(task, &_randomVals);

	//Allocate memory for the track after being randomly shifted
	_randomlyShiftedTrack.AllocateEmptyTrack(_nSteps);

	//Allocate memory for the number of steps within the bounding box
	cudaMallocManaged(&_numInVoxel,sizeof(int)); 

	//Allocate memory to store the StepIDs of the steps within spheres
	cudaMalloc(&_inSphereTrackId,_nSteps*sizeof(int));
}

void VoxelConstrainedSphereMethod::GenerateRandomXYShift(const ThreadTask &task, float** randomVals)
{
	cudaMalloc(randomVals,2*sizeof(float)*task.GetNOversamples()); 

	if (_randomShifts == true)
	{
		//Create the random generator
		curandGenerator_t randGenerator;
		curandCreateGenerator(&randGenerator,CURAND_RNG_PSEUDO_DEFAULT);

		//Seed the generator
		curandSetPseudoRandomGeneratorSeed(randGenerator,task.GetRandomSeed());

		//Make random numbers, and then destroy the generator
		curandGenerateUniform(randGenerator,*randomVals,2*task.GetNOversamples());
		curandDestroyGenerator(randGenerator);
	}
	else //set the random vals to 0.5 (no shift) if not requested
	{
		VoxelConstrainedSphereMethodKernel::SetRandomValsToHalf<<<_suggestedCudaBlocks,_suggestedCudaThreads>>>(*randomVals, 2*task.GetNOversamples());
	}
	cudaDeviceSynchronize();
}

//Called repeatedly for each track oversample
void VoxelConstrainedSphereMethod::ProcessTrack(Track track, VolumeEdepPair& edepsInTarget)
{ 
	//New track. Zero values
	SimulationMethodKernel::ZeroInt<<<1,1>>>(_numInVoxel);
	SimulationMethodKernel::ZeroInt<<<1,1>>>(edepsInTarget.numElements);

	//Filter and score the tracks
		//Filter the tracks that are in the scoring voxel (box)
	VoxelConstrainedSphereMethodKernel::FilterInScoringBox<<<_suggestedCudaBlocks,_suggestedCudaThreads>>>(_sphericalGeometry,_randomVals,track,_randomlyShiftedTrack,_nSteps,_numInVoxel,_oversampleIterationNumber);	
		//Filter the tracks that land inside of a scoring sphere
	VoxelConstrainedSphereMethodKernel::FilterTrackInSphere<<<_suggestedCudaBlocks,_suggestedCudaThreads>>>(_sphericalGeometry,_randomlyShiftedTrack,_numInVoxel,edepsInTarget.numElements,_inSphereTrackId); 
		//Score the the tracks which reside in a sphere
	VoxelConstrainedSphereMethodKernel::ScoreTrackInSphere<<<_suggestedCudaBlocks,_suggestedCudaThreads>>>(_sphericalGeometry,_randomlyShiftedTrack,edepsInTarget.numElements,_inSphereTrackId,edepsInTarget);

	_oversampleIterationNumber++;
}

//Called at the end of processing a track
void VoxelConstrainedSphereMethod::FreeTrackProcess()
{ 
	//Free directly allocated memory
	cudaFree(_inSphereTrackId);
	cudaFree(_randomVals);
	cudaFree(_numInVoxel);

	//Free my classes
	_randomlyShiftedTrack.Free();
}

//Called at the end of processing all of the tracks
void VoxelConstrainedSphereMethod::Free()
{ 

}

//
//Kernel definitions
//

__global__ void VoxelConstrainedSphereMethodKernel::SetRandomValsToHalf(float* randomVals, int n)
{
	//Determine index and strid
    int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;

	//Loop over all random vals and set to 0.5
	for (int i = index; i < n; i+=stride)
	{
		randomVals[i] = 0.5;
	}
}

__global__ void VoxelConstrainedSphereMethodKernel::FilterInScoringBox(SphericalGeometry geometry, float* randomVals, Track inputTrack, Track outputTrack, int numElements, int *numElementsCompacted, int oversampleIterationNumber)
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
	
	//Shift a random float from 0-1 into the range from -halflength to +halflength
	double x_shift = ((randomVals[(oversampleIterationNumber*2)]*geometry.scoringRegionLength)-geometry.scoringRegionHalfLength);
	double y_shift = ((randomVals[(oversampleIterationNumber*2+1)]*geometry.scoringRegionLength)-geometry.scoringRegionHalfLength);

	//Loop over all the energy deposition points
	for (int i = index; i < numElements; i+=stride)
	{
		//Apply the random shift
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
		//we don't have to check Z, the tracks are generated so they are never outside in Z
		if (abs(x_shifted) < geometry.scoringRegionHalfLength  && abs(y_shifted) < geometry.scoringRegionHalfLength) 
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

			//Copy the track inside the box over to the new array
			outputTrack.x[outputIndex] = x_shifted;
			outputTrack.y[outputIndex] = y_shifted;
			outputTrack.z[outputIndex] = inputTrack.z[i];
			outputTrack.edep[outputIndex] = inputTrack.edep[i];
		}
		__syncthreads();
	}
}

__global__ void VoxelConstrainedSphereMethodKernel::FilterTrackInSphere(SphericalGeometry geometry, Track inputTrack, int *numElements, int *numElementsCompacted, int *trackIdInSphere)
{
	//move all of the variable definitions out of the for loop
	double distFromNearestSphereX, distFromNearestSphereY, distFromNearestSphereZ, dist;

	//Determine index and stride
 	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;

	//Counters for shared memory atomics
	int localPosition;
	__shared__ int localIndexCounter;

	//Pre-calculate values
	double sphereRadiusMag = geometry.sphereRadius*geometry.sphereRadius; 

	//Loop over all the energy deposition points
	for (long i = index; i < *numElements; i+=stride)
	{		
		//Take the position relative to the greatest sphere offset. Divide by the number of sphere diameters away, and round to nearest int to find the index in each axis
		//Then multiple the index by the sphere diameter, and subtract the current position to get the distance away

		double distFromFurthestSphereX = inputTrack.x[i]-geometry.greatestSphereOffset;
		double distFromFurthestSphereY = inputTrack.y[i]-geometry.greatestSphereOffset;
		double distFromFurthestSphereZ = inputTrack.z[i]-geometry.greatestSphereOffset;

		distFromNearestSphereX = ((llrint(distFromFurthestSphereX/geometry.sphereDiameter))*geometry.sphereDiameter)-distFromFurthestSphereX;
		distFromNearestSphereY = ((llrint(distFromFurthestSphereY/geometry.sphereDiameter))*geometry.sphereDiameter)-distFromFurthestSphereY;
		distFromNearestSphereZ = ((llrint(distFromFurthestSphereZ/geometry.sphereDiameter))*geometry.sphereDiameter)-distFromFurthestSphereZ; 

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

__global__ void VoxelConstrainedSphereMethodKernel::ScoreTrackInSphere(SphericalGeometry geometry, Track inputTrack, int *numElements, int *trackIdInSphere, VolumeEdepPair outputPair)
{
	//move all of the variable definitions out of the for loop
	uint64_t xIndex, yIndex, zIndex, sphereHitIndex;

	//Pre-calculate some values
	double sphereDiameter = geometry.sphereDiameter;
	double linealDenominator = (2./3.)*sphereDiameter; 

	//Determine index and stride
 	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;

	//Loop over all the energy deposition points
	for (uint64_t i = index; i < *numElements; i+=stride)
	{	
		//Take the position relative to the greatest sphere offset. Divide by the number of sphere diameters away,
		//and round to nearest int to find the index in each axis
		xIndex = llrint((inputTrack.x[trackIdInSphere[i]]-geometry.greatestSphereOffset)/geometry.sphereDiameter);
		yIndex = llrint((inputTrack.y[trackIdInSphere[i]]-geometry.greatestSphereOffset)/geometry.sphereDiameter);
		zIndex = llrint((inputTrack.z[trackIdInSphere[i]]-geometry.greatestSphereOffset)/geometry.sphereDiameter);

		//Determine the Index of the sphere hit
		sphereHitIndex = xIndex + yIndex*geometry.numSpheresLinear+ zIndex*geometry.numSpheresLinear*geometry.numSpheresLinear; //Keep in mind that for the index it starts counting at zero

		//Write to volumeID and edepOutput
		outputPair.volume[i] = sphereHitIndex;
		outputPair.edep[i] = inputTrack.edep[trackIdInSphere[i]]/linealDenominator; //this should be ev/nm which is same a kev/um
	}
}

