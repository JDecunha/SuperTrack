#include "VolumeEdepPair.cuh"
#include <iostream>

VolumeEdepPair::VolumeEdepPair() {}

void VolumeEdepPair::Allocate(uint64_t numInputElements)
{
	cudaMalloc(&volume,numInputElements*sizeof(uint64_t));
	cudaMalloc(&edep,numInputElements*sizeof(double));
	cudaMallocManaged(&numElements,sizeof(int));

	//Set number of elements from number of elements allocated
	*numElements = numInputElements;
	//Zero the volume edep pairs
	cudaMemset(volume,0,numInputElements*sizeof(uint64_t));
	cudaMemset(edep,0,numInputElements*sizeof(double));
}

void VolumeEdepPair::Free()
{
	cudaFree(volume);
	cudaFree(edep);
	cudaFree(numElements);
}



