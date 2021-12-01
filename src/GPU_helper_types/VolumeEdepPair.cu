#include "VolumeEdepPair.cuh"
#include <iostream>

VolumeEdepPair::VolumeEdepPair() {}

void VolumeEdepPair::Allocate(uint64_t numInputElements)
{
	cudaMalloc(&volume,numInputElements*sizeof(uint64_t));
	cudaMalloc(&edep,numInputElements*sizeof(double));
	cudaMallocManaged(&numElements,sizeof(int));

	*numElements = numInputElements;
}

void VolumeEdepPair::Free()
{
	cudaFree(volume);
	cudaFree(edep);
	cudaFree(numElements);
}



