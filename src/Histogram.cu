#include "Histogram.cuh"
#include "utils.hh"

//
//Histogram definitions
//

Histogram::Histogram(int nbins, float binLower, float binUpper,std::string type="log")
{
	_nbins = nbins;
	_binLower = binLower;
	_binUpper = binUpper;
	_type = type;
}

void Histogram::GenerateLogHistogram()
{
	//Get the device Id for active GPU
	int deviceId;	cudaGetDevice(&deviceId);   

	//Fill the log bins and send to the device
	cudaMallocManaged(&_binEdges, (_nbins+1)*sizeof(double));
	LogSpace(_binLower,_binUpper,_nbins,_binEdges);
	cudaMemPrefetchAsync(_binEdges,(_nbins+1)*sizeof(double),deviceId);

	//TODO: Change to unmanaged memory later
	cudaMallocManaged(&_histogramVals,_nbins*sizeof(int));
	cudaMallocManaged(&_histogramValsAccumulated,_nbins*sizeof(int));

	//Set arrays to zero
	cudaMemset(_histogramVals,0,_nbins*sizeof(int));
	cudaMemset(_histogramValsAccumulated,0,_nbins*sizeof(int));
}

void Histogram::Allocate(VolumeEdepPair targetEdeps)
{
	if (_type == "log")
	{
		GenerateLogHistogram();
	}
	else
	{
		std::cout << "Non-log histogram not yet supported in SuperTrack" << std::endl;
	}

	//Allocate VolumeID-Edep pairs
	sortedEdeps.Allocate(*(targetEdeps.numElements));
	reducedEdeps.Allocate(*(targetEdeps.numElements));

	//Allocate my CubStorageBuffers
	sortBuffer = CubStorageBuffer::AllocateCubSortBuffer(targetEdeps);
	reduceBuffer = CubStorageBuffer::AllocateCubReduceBuffer(targetEdeps);
	histogramBuffer = CubStorageBuffer::AllocateCubHistogramBuffer(targetEdeps,_histogramVals,_binEdges,_nbins);

}

void Histogram::Free()
{
	cudaFree(_binEdges);
	cudaFree(_histogramVals);
	cudaFree(_histogramValsAccumulated);

	sortedEdeps.Free();
	reducedEdeps.Free();

	sortBuffer.Free();
	reduceBuffer.Free();
	histogramBuffer.Free();
}

void Histogram::Accumulate()
{
	HistogramKernel::AccumulateHistogramVals<<<4,32>>>(_histogramVals,_histogramValsAccumulated,_nbins);
}

void Histogram::SortReduceAndAddToHistogram(VolumeEdepPair targetEdeps)
{
	HistogramKernel::SortReduceAndAddToHistogramKernel<<<1,1>>>(sortBuffer,reduceBuffer,histogramBuffer,targetEdeps,sortedEdeps,reducedEdeps, _nbins,_histogramVals, _binEdges, reductionOperator);
	Accumulate();
}

//
//HistogramKernel definitions
//

__global__ void HistogramKernel::AccumulateHistogramVals(int* temp, int* accumulated,int N)
{
	//Determine index and stride
 	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;

	for (int i = index; i < N; i+=stride)
	{
		accumulated[i] = accumulated[i]+temp[i];
	}
}

__global__ void HistogramKernel::SortReduceAndAddToHistogramKernel(CubStorageBuffer sortBuffer, CubStorageBuffer reduceBuffer, CubStorageBuffer histogramBuffer, VolumeEdepPair edepsInTarget, VolumeEdepPair sortedEdeps, VolumeEdepPair reducedEdeps, int nbins,int* histogramVals, double* logBins, CUBAddOperator reductionOperator)
{
	//Sort the edep volume pairs
	cub::DeviceRadixSort::SortPairs(sortBuffer.storage,sortBuffer.size,edepsInTarget.volume,sortedEdeps.volume,edepsInTarget.edep,sortedEdeps.edep,*(edepsInTarget.numElements));
	// reduce the energy depositions
	cub::DeviceReduce::ReduceByKey(reduceBuffer.storage,reduceBuffer.size, sortedEdeps.volume, reducedEdeps.volume, sortedEdeps.edep, reducedEdeps.edep, reducedEdeps.numElements, reductionOperator, *(edepsInTarget.numElements));
	//Create the histogram
	cub::DeviceHistogram::HistogramRange(histogramBuffer.storage,histogramBuffer.size,reducedEdeps.edep,histogramVals,nbins+1,logBins,*reducedEdeps.numElements);
}



