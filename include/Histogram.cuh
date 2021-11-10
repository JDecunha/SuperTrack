#pragma once

#include <string>
#include "CubStorageBuffer.cuh"
#include "VolumeEdepPair.cuh"

class Histogram
{
	public:
		Histogram(int nbins, float binLower, float binUpper, std::string type);

		void Accumulate();
		void Allocate();
		void Free();
		
		int _nbins;
		float _binLower, _binUpper;
		double *_binEdges;
		int *_histogramVals, *_histogramValsAccumulated;

	private:
		void GenerateLogHistogram();

		VolumeEdepPair targetEdeps, sortedEdeps, reducedEdeps;
		CubStorageBuffer sortBuffer, reduceBuffer, histogramBuffer;

		std::string _type;
};

//CUDA Limitation: Kernels cannot belong to classes
//So to keep things organized I'm putting them into a class associated namespace
namespace HistogramKernel
{
	__global__ void AccumulateHistogramVals(int* temp, int* accumulated,int N);
};