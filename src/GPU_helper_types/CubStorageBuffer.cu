#include "CubStorageBuffer.cuh"
#include <cub/cub.cuh>

CubStorageBuffer::CubStorageBuffer()
{
	storage = NULL;
	size = 0;
}

void CubStorageBuffer::Free()
{
	cudaFree(storage);
}

CubStorageBuffer CubStorageBuffer::AllocateCubSortBuffer(VolumeEdepPair edepPairList)
{
	//Create the buffer with default constructor
	CubStorageBuffer returnBuffer = CubStorageBuffer();

	//Call the CUB function to determine memory constraints, then malloc
	cub::DeviceRadixSort::SortPairs(returnBuffer.storage,returnBuffer.size,edepPairList.volume,edepPairList.volume,edepPairList.edep,edepPairList.edep,*(edepPairList.numElements));
	cudaMalloc(&returnBuffer.storage,returnBuffer.size); 

	return returnBuffer;
}

CubStorageBuffer CubStorageBuffer::AllocateCubReduceBuffer(VolumeEdepPair edepPairList)
{
	//Create the buffer with default constructor
	CubStorageBuffer returnBuffer = CubStorageBuffer();
	//Generic Reduction Operator
	CUBAddOperator reductionOperator;

	//Call the CUB function to determine memory constraints, then malloc
	cub::DeviceReduce::ReduceByKey(returnBuffer.storage,returnBuffer.size, edepPairList.volume, edepPairList.volume, edepPairList.edep, edepPairList.edep, edepPairList.numElements, reductionOperator, *(edepPairList.numElements));
	cudaMalloc(&returnBuffer.storage,returnBuffer.size); 

	return returnBuffer;
}

CubStorageBuffer CubStorageBuffer::AllocateCubHistogramBuffer(VolumeEdepPair edepPairList, int* histogramVals, double* logBins, int nbins)
{
	//Create the buffer with default constructor
	CubStorageBuffer returnBuffer = CubStorageBuffer();

	//Call the CUB function to determine memory constraints, then malloc
	cub::DeviceHistogram::HistogramRange(returnBuffer.storage,returnBuffer.size, edepPairList.edep,histogramVals,nbins+1,logBins,*(edepPairList.numElements));
	cudaMalloc(&returnBuffer.storage,returnBuffer.size); 

	return returnBuffer;
}