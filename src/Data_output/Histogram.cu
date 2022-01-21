#include "Histogram.cuh"
#include "utils.hh"

//
// Histogram definitions
//

//
// Constructor and destructor related functions
//

Histogram::Histogram(int nbins, float binLower, float binUpper,std::string type="log")
{
	_nbins = nbins;
	_binLower = binLower;
	_binUpper = binUpper;
	_type = type;

	//Values have been set now allocate memory
	Allocate();
}

Histogram::Histogram(const INIReader& reader)
{
	_type = reader.Get("Histogram","Type","");
	_binLower = reader.GetFloat("Histogram","BinLower",0);
	_binUpper = reader.GetFloat("Histogram","BinUpper",0);
	_nbins = reader.GetInteger("Histogram","NBins",0);

	//Values have been set now allocate memory
	Allocate();

	//Additional values not required to allocate the GPU histogram
	_suggestedCudaAccumulateBlocks = reader.GetReal("Histogram","SuggestedCudaAccumulateBlocks",4);
	_suggestedCudaAccumulateThreads = reader.GetReal("Histogram","SuggestedCudaAccumulateThreads",32);
}

Histogram::~Histogram()
{
	//Free all the memory from when the object was initially constructed
	cudaFree(_binEdges);
	cudaFree(_histogramVals);
	cudaFree(_histogramValsAccumulated);
}

void Histogram::Allocate()
{
	if (_type == "log")
	{
		GenerateLogHistogram();
		InitializeCPULogHistogram(); //Generate the CPU histogram
	}
	else
	{
		std::cout << "Non-log histogram not yet supported in SuperTrack" << std::endl; abort();
	}
}

void Histogram::GenerateLogHistogram()
{
	//Get the device Id for active GPU
	int deviceId;	cudaGetDevice(&deviceId);   

	//Fill the log bins and send to the device
	cudaMallocManaged(&_binEdges, (_nbins+1)*sizeof(double));
	utils::LogSpace(_binLower,_binUpper,_nbins,_binEdges);
	cudaMemPrefetchAsync(_binEdges,(_nbins+1)*sizeof(double),deviceId);

	//TODO: Change to unmanaged memory later
	cudaMallocManaged(&_histogramVals,_nbins*sizeof(int));
	cudaMallocManaged(&_histogramValsAccumulated,_nbins*sizeof(int));

	//Set arrays to zero
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

void Histogram::Accumulate()
{
	HistogramKernel::AccumulateHistogramVals<<<_suggestedCudaAccumulateBlocks,_suggestedCudaAccumulateThreads>>>(_histogramVals,_histogramValsAccumulated,_nbins);
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
	//Reduce the energy depositions
	cub::DeviceReduce::ReduceByKey(reduceBuffer.storage,reduceBuffer.size, sortedEdeps.volume, reducedEdeps.volume, sortedEdeps.edep, reducedEdeps.edep, reducedEdeps.numElements, reductionOperator, *(edepsInTarget.numElements));
	//Create the histogram
	cub::DeviceHistogram::HistogramRange(histogramBuffer.storage,histogramBuffer.size,reducedEdeps.edep,histogramVals,nbins+1,logBins,*reducedEdeps.numElements);
}



//
//CPU Histogram related
//

void Histogram::InitializeCPULogHistogram()
{
	//Initialize the histogram with appropriate values
	 _CPU_histogram = TH1D("Lineal energy histogram", "y*f(y)", _nbins, _binLower,_binUpper); 
	//Make the bins logarithmic
	CPUHistogramUtils::BinLogX(&_CPU_histogram);
}

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



