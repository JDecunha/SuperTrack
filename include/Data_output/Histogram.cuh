#pragma once

//INIH
#include "INIReader.h"
//SuperTrack
#include "CubAddOperator.cuh"
#include "CubStorageBuffer.cuh"
#include "VolumeEdepPair.cuh"
//ROOT
#include "TH1D.h"
//STD
#include <string>


class Histogram
{
	public:
		//Can construct directly or with user inputs
		Histogram(int nbins, float binLower, float binUpper, std::string type, int suggestedCudaAccumulateBlocks = 4, int suggestedCudaAccumulatedThreads = 32);
		Histogram(const INIReader& reader);

		//Custom destructor to free memory allocations
		~Histogram();

		//Object should never be assigned or copied
		Histogram(Histogram const&) = delete;
		void operator=(Histogram const&) = delete;

		//These get called when the track is being analyzed so only values that change with each track are allocated and freed
		void AllocateTrackProcess(VolumeEdepPair targetEdeps);
		void FreeTrackProcess();

		//Take the edeps and sort, reduce, and add them to the histogram
		void SortReduceAndAddToHistogram(VolumeEdepPair targetEdeps);

		//Outputting the histogram
		void Print();

		//Return the CPU histogram
		TH1D GetCPUHistogram() { return _CPU_histogram; }

		//Helper functions
		static TH1D ReduceVector(std::vector<TH1D>& toReduce);

	private:
		void ConstructBins();

		//Functions for the GPU histogram
		void ConstructGPUHistogram();
		void AccumulateGPUHistogram();
		
		//Functions for the CPU histogram
		TH1D _CPU_histogram;
		void ConstructCPUHistogram();
		void AccumulateCPUHistogram();	//Add the current values from the GPU histogram to the CPU histogram

		//Attributes
		std::string _type;
		int _nbins;
		float _binLower, _binUpper;
		int _suggestedCudaAccumulateThreads;
		int _suggestedCudaAccumulateBlocks;

		//Attributes in CudaMalloced memory
		double *_binEdges;
		int *_histogramVals, *_histogramValsAccumulated;

		//Helper classes and storage structs
		VolumeEdepPair sortedEdeps, reducedEdeps;
		CubStorageBuffer sortBuffer, reduceBuffer, histogramBuffer;
		CUBAddOperator reductionOperator;
};

//CUDA Limitation: Kernels cannot belong to classes
//So to keep things organized I'm putting them into a class associated namespace
namespace HistogramKernel
{
	__global__ void AccumulateHistogramVals(int* temp, int* accumulated,int N);
	__global__ void SortReduceAndAddToHistogramKernel(CubStorageBuffer sortBuffer, CubStorageBuffer reduceBuffer, CubStorageBuffer histogramBuffer, VolumeEdepPair edepsInTarget, VolumeEdepPair sortedEdeps, VolumeEdepPair reducedEdeps, int nbins, int* histogramVals, double* logBins, CUBAddOperator reductionOperator);
};