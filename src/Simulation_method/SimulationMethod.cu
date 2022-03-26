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

//
// Kernel definitions
//

__global__ void SimulationMethodKernel::ZeroInt(int* toZero)
{
	*toZero = 0;
}


