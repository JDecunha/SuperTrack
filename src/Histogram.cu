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

void Histogram::Allocate()
{
	if (_type == "log")
	{
		GenerateLogHistogram();
	}
	else
	{
		std::cout << "Non-log histogram not yet supported in SuperTrack" << std::endl;
	}
}

void Histogram::Free()
{
	cudaFree(_binEdges);
	cudaFree(_histogramVals);
	cudaFree(_histogramValsAccumulated);
}


void Histogram::GenerateLogHistogram()
{
	//Get the device Id for active GPU
	int deviceId;
	cudaGetDevice(&deviceId);   

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

void Histogram::Accumulate()
{
	HistogramKernel::AccumulateHistogramVals<<<4,32>>>(_histogramVals,_histogramValsAccumulated,_nbins);
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



