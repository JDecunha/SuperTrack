//SuperTrack
#include "SimulationMethod.hh"
#include "Track.cuh"
#include "VolumeEdepPair.cuh"
//CUDA Libraries
#include <cuda.h>
#include <curand.h>


SimulationMethod::SimulationMethod(const INIReader& macroReader)
{
	_macroReader = macroReader;
}

void SimulationMethod::GenerateRandomXYShift(const ThreadTask &task, float** randomVals)
{
	cudaMalloc(randomVals,2*sizeof(float)*task.GetNOversamples()); 
	
	//Random number generation on GPU
	curandGenerator_t randGenerator;
	curandCreateGenerator(&randGenerator,CURAND_RNG_PSEUDO_DEFAULT);
	//TODO: Change this to collate the numbers rather than add the threadID to randomSeed
	curandSetPseudoRandomGeneratorSeed(randGenerator,task.GetRandomSeed()+task.GetThreadID());
	curandGenerateUniform(randGenerator,*randomVals,2*task.GetNOversamples());
	curandDestroyGenerator(randGenerator);
	cudaDeviceSynchronize();
}

//
// Kernel definitions
//

__global__ void SimulationMethodKernel::ZeroInt(int* toZero)
{
	*toZero = 0;
}


