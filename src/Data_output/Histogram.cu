#include "Histogram.cuh"
#include "utils.hh"

//
// Histogram definitions
//

//References:
//1.) CUB Device Histogram reference page: https://nvlabs.github.io/cub/structcub_1_1_device_histogram.html#a11ec6d941fb6779c2a4d124b6f5b0813
//2.) TH1 reference page: https://root.cern.ch/doc/master/classTH1.html#a82f16fb9b9a11c97f1c9177e9e996fc7

//
// Constructor and destructor related functions
//

Histogram::Histogram(int nbins, float binLower, float binUpper,std::string type="log", int suggestedCudaAccumulateBlocks, int suggestedCudaAccumulatedThreads)
{
	_nbins = nbins;
	_binLower = binLower;
	_binUpper = binUpper;
	_type = type;
	_suggestedCudaAccumulateBlocks = suggestedCudaAccumulateBlocks;
	_suggestedCudaAccumulateThreads = suggestedCudaAccumulatedThreads;

	ConstructBins(); //Creates the bin edges used both by the CPU and GPU histogram
	ConstructCPUHistogram();
	ConstructGPUHistogram();
}

Histogram::Histogram(const INIReader& reader)
{
	//Read in values from macro file
	_type = reader.Get("Histogram","Type","");
	_binLower = reader.GetFloat("Histogram","BinLower",0);
	_binUpper = reader.GetFloat("Histogram","BinUpper",0);
	_nbins = reader.GetInteger("Histogram","NBins",0);
	_suggestedCudaAccumulateBlocks = reader.GetReal("Histogram","SuggestedCudaAccumulateBlocks",4);
	_suggestedCudaAccumulateThreads = reader.GetReal("Histogram","SuggestedCudaAccumulateThreads",32);

	ConstructBins(); //Creates the bin edges used both by the CPU and GPU histogram
	ConstructCPUHistogram();
	ConstructGPUHistogram();
}

Histogram::~Histogram()
{
	//Free all the memory from when the object was initially constructed
	cudaFree(_binEdges);
	cudaFree(_histogramVals);
	cudaFree(_histogramValsAccumulated);
}

void Histogram::ConstructBins()
{
	//Allocate memory for the bin edges
	cudaMallocManaged(&_binEdges, (_nbins+1)*sizeof(double));

	//Set the bin edges based on histogram type
	if (_type == "log")
	{
		utils::LogSpace(_binLower,_binUpper,_nbins,_binEdges);
	}
	else if (_type == "lin")
	{
		utils::LinSpace(_binLower,_binUpper,_nbins,_binEdges);
	}
	else
	{
		std::cout << "Unknown histogram type: " << _type << " available types are lin and log." << std::endl; abort();
	}
	
}

void Histogram::ConstructCPUHistogram()
{
	//Initialize the histogram with appropriate values
	_CPU_histogram = TH1D("Lineal energy histogram", "N(y)", _nbins, _binEdges); 
}

void Histogram::ConstructGPUHistogram()
{
	//Get the device Id for active GPU
	int deviceId;	cudaGetDevice(&deviceId);   

	//Pre-fetch the bin edges to the GPU
	cudaMemPrefetchAsync(_binEdges,(_nbins+1)*sizeof(double),deviceId);

	//Create the histogram bin values and accumulated bin values
	//TODO: Change to unmanaged memory later
	cudaMallocManaged(&_histogramVals,_nbins*sizeof(int));
	cudaMallocManaged(&_histogramValsAccumulated,_nbins*sizeof(int));

	//Zero bin values
	cudaMemset(_histogramVals,0,_nbins*sizeof(int));
	cudaMemset(_histogramValsAccumulated,0,_nbins*sizeof(int));
}

//
// Functions called during track processing
//

void Histogram::AllocateTrackProcess(VolumeEdepPair targetEdeps)
{
	//Allocate VolumeID-Edep pairs
	sortedEdeps.Allocate(*(targetEdeps.numElements));
	reducedEdeps.Allocate(*(targetEdeps.numElements));

	//Allocate my CubStorageBuffers
	sortBuffer = CubStorageBuffer::AllocateCubSortBuffer(targetEdeps);
	reduceBuffer = CubStorageBuffer::AllocateCubReduceBuffer(targetEdeps);
	histogramBuffer = CubStorageBuffer::AllocateCubHistogramBuffer(targetEdeps,_histogramVals,_binEdges,_nbins);
}

void Histogram::SortReduceAndAddToHistogram(VolumeEdepPair targetEdeps)
{
	HistogramKernel::SortReduceAndAddToHistogramKernel<<<1,1>>>(sortBuffer,reduceBuffer,histogramBuffer,targetEdeps,sortedEdeps,reducedEdeps, _nbins,_histogramVals, _binEdges, reductionOperator);
	AccumulateGPUHistogram();
}

void Histogram::AccumulateGPUHistogram()
{
	HistogramKernel::AccumulateHistogramVals<<<_suggestedCudaAccumulateBlocks,_suggestedCudaAccumulateThreads>>>(_histogramVals,_histogramValsAccumulated,_nbins);
}

void Histogram::FreeTrackProcess()
{
	//End of track process, accumulate CPU histogram
	AccumulateCPUHistogram();

	//Then free memory
	sortedEdeps.Free();
	reducedEdeps.Free();

	sortBuffer.Free();
	reduceBuffer.Free();
	histogramBuffer.Free();

	//Don't need to reallocate these values, just zero them for the next track
	cudaMemset(_histogramVals,0,_nbins*sizeof(int));
	cudaMemset(_histogramValsAccumulated,0,_nbins*sizeof(int));
}

//
//HistogramKernel definitions
//

__global__ void HistogramKernel::SortReduceAndAddToHistogramKernel(CubStorageBuffer sortBuffer, CubStorageBuffer reduceBuffer, CubStorageBuffer histogramBuffer, VolumeEdepPair edepsInTarget, VolumeEdepPair sortedEdeps, VolumeEdepPair reducedEdeps, int nbins,int* histogramVals, double* logBins, CUBAddOperator reductionOperator)
{
	//Sort the edep volume pairs
	cub::DeviceRadixSort::SortPairs(sortBuffer.storage,sortBuffer.size,edepsInTarget.volume,sortedEdeps.volume,edepsInTarget.edep,sortedEdeps.edep,*(edepsInTarget.numElements));

	cudaDeviceSynchronize();

	//Reduce the energy depositions
	cub::DeviceReduce::ReduceByKey(reduceBuffer.storage,reduceBuffer.size, sortedEdeps.volume, reducedEdeps.volume, sortedEdeps.edep, reducedEdeps.edep, reducedEdeps.numElements, reductionOperator, *(edepsInTarget.numElements));

	cudaDeviceSynchronize();

	//Create the histogram
	cub::DeviceHistogram::HistogramRange(histogramBuffer.storage,histogramBuffer.size,reducedEdeps.edep,histogramVals,nbins+1,logBins,*reducedEdeps.numElements);

	cudaDeviceSynchronize();
}

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

//
//CPU Histogram related
//

void Histogram::AccumulateCPUHistogram()
{
	//Read out GPU histogram and fill CPU histogram
	for (int i = 0; i < _nbins; i++)
	{
		//Fill each bin with appropriate value
		//i+1 because 0th bin is the underflow value
		_CPU_histogram.AddBinContent(i+1,_histogramValsAccumulated[i]);
	}
}

TH1D Histogram::ReduceVector(std::vector<TH1D>& toReduce)
{
	TH1D returnVal = toReduce.back();
	toReduce.pop_back(); //Remove the last element because we just took it

	for (TH1D histogram : toReduce) //Add the other histograms to this one
	{
		returnVal.Add(&histogram);
	}

	return returnVal;
}

//
// Misc
//

void Histogram::Print()
{
	int number_of_values_in_histogram = 0;

	//TODO: place a cudamemcpy here, so that the memory transfer is faster
				
	//Read out histogram
	for (int i = 0; i < _nbins; i++)
	{
		number_of_values_in_histogram += _histogramValsAccumulated[i];
		std::cout << "Bin: " << _binEdges[i] << " Counts: " << _histogramValsAccumulated[i] << std::endl;
	}
	std::cout << number_of_values_in_histogram << std::endl;
}



