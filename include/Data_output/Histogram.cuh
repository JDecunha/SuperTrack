#pragma once

#include <string>

//CUB
#include "CubAddOperator.cuh"
#include "CubStorageBuffer.cuh"
#include "VolumeEdepPair.cuh"


class Histogram
{
	public:
		Histogram(int nbins, float binLower, float binUpper, std::string type);

		void Accumulate();
		void Allocate(VolumeEdepPair targetEdeps);
		void Free();
		void SortReduceAndAddToHistogram(VolumeEdepPair targetEdeps);
		
		int _nbins;
		float _binLower, _binUpper;
		double *_binEdges;
		int *_histogramVals, *_histogramValsAccumulated;

		VolumeEdepPair sortedEdeps, reducedEdeps;
		CubStorageBuffer sortBuffer, reduceBuffer, histogramBuffer;
		CUBAddOperator reductionOperator;

	private:
		void GenerateLogHistogram();

		

		std::string _type;
};

//CUDA Limitation: Kernels cannot belong to classes
//So to keep things organized I'm putting them into a class associated namespace
namespace HistogramKernel
{
	__global__ void AccumulateHistogramVals(int* temp, int* accumulated,int N);
	__global__ void SortReduceAndAddToHistogramKernel(CubStorageBuffer sortBuffer, CubStorageBuffer reduceBuffer, CubStorageBuffer histogramBuffer, VolumeEdepPair edepsInTarget, VolumeEdepPair sortedEdeps, VolumeEdepPair reducedEdeps, int nbins,int* histogramVals, double* logBins, CUBAddOperator reductionOperator);
};