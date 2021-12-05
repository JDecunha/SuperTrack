//SuperTrack
#include "SimulationMethod.hh"
#include "Track.cuh"
#include "VolumeEdepPair.cuh"
//CUDA Libraries
#include <cuda.h>
#include <curand.h>
#include <iostream>


SimulationMethod::SimulationMethod(const INIReader& macroReader)
{
	_macroReader = macroReader;
}

void SimulationMethod::GenerateRandomXYShift(const ThreadTask &task, float** randomVals)
{
	cudaMalloc(randomVals,2*sizeof(float)*task.GetNOversamples()); 
	
	//Create the random generator
	curandGenerator_t randGenerator;
	curandCreateGenerator(&randGenerator,CURAND_RNG_PSEUDO_DEFAULT);

	//Seed the generator
	curandSetPseudoRandomGeneratorSeed(randGenerator,task.GetRandomSeed());

	//Make random numbers, and then destroy the generator
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


