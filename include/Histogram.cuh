#pragma once

#include <string>

class Histogram
{
	public:
		Histogram(int nbins, float binLower, float binUpper, std::string type);
		~Histogram();
		void Accumulate();
		
		int _nbins;
		float _binLower, _binUpper;
		double *_binEdges;
		int *_histogramVals, *_histogramValsAccumulated;

	private:
		void GenerateLogHistogram();
};

//CUDA Limitation: Kernels cannot belong to classes
//So to keep things organized I'm putting them into a class associated namespace
namespace HistogramKernel
{
	__global__ void AccumulateHistogramVals(int* temp, int* accumulated,int N);
};