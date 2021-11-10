#pragma once

#include "CubStorageBuffer.cuh"
#include "SuperTrackTypes.cuh"
#include "VolumeEdepPair.cuh"

class CubStorageBuffer
{
	public:
		CubStorageBuffer();
		~CubStorageBuffer();

		static CubStorageBuffer AllocateCubSortBuffer(VolumeEdepPair edepPairList, uint64_t nVals);
		static CubStorageBuffer AllocateCubReduceBuffer(VolumeEdepPair edepPairList, uint64_t nVals);
		static CubStorageBuffer AllocateCubHistogramBuffer(VolumeEdepPair edepPairList, uint64_t nVals, int* histogramVals, double* logBins, int nbins);
		
		void* storage;
		size_t size;
};