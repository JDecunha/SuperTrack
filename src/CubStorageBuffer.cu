#include "CubStorageBuffer.cuh"
#include <cub/cub.cuh>

CubStorageBuffer::CubStorageBuffer()
{
	storage = NULL;
	size = 0;
}

CubStorageBuffer::~CubStorageBuffer()
{
	cudaFree(&storage);
	cudaFree(&size);
}

CubStorageBuffer CubStorageBuffer::AllocateCubSortBuffer(VolumeEdepPair edepPairList, uint64_t nVals)
{
	//Create the buffer with default constructor
	CubStorageBuffer returnBuffer = CubStorageBuffer();

	//Call the CUB function to determine memory constraints, then malloc
	cub::DeviceRadixSort::SortPairs(returnBuffer.storage,returnBuffer.size,edepPairList.volume,edepPairList.volume,edepPairList.edep,edepPairList.edep,nVals);
	cudaMalloc(&returnBuffer.storage,returnBuffer.size); 

	return returnBuffer;
}

CubStorageBuffer CubStorageBuffer::AllocateCubReduceBuffer(VolumeEdepPair edepPairList, uint64_t nVals)
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

CubStorageBuffer CubStorageBuffer::AllocateCubHistogramBuffer(VolumeEdepPair edepPairList, uint64_t nVals, int* histogramVals, double* logBins, int nbins)
{
	//Create the buffer with default constructor
	CubStorageBuffer returnBuffer = CubStorageBuffer();

	//Call the CUB function to determine memory constraints, then malloc
	cub::DeviceHistogram::HistogramRange(returnBuffer.storage,returnBuffer.size, edepPairList.edep,histogramVals,nbins+1,logBins,nVals);
	cudaMalloc(&returnBuffer.storage,returnBuffer.size); 

	return returnBuffer;
}