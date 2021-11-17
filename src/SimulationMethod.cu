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

//
// Kernel definitions
//

__global__ void SimulationMethodKernel::ZeroInt(int* toZero)
{
	*toZero = 0;
}

void SimulationMethodKernel::GenerateRandomXYShift(const std::tuple<Int_t,Int_t,Int_t,TString> &input, float **randomVals, const int &nSamples, const long &random_seed)
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