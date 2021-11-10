#pragma once

#include "CubStorageBuffer.cuh"
#include "SuperTrackTypes.cuh"
#include "VolumeEdepPair.cuh"

class CubStorageBuffer
{
	public:
		CubStorageBuffer();
		void Free();

		static CubStorageBuffer AllocateCubSortBuffer(VolumeEdepPair edepPairList);
		static CubStorageBuffer AllocateCubReduceBuffer(VolumeEdepPair edepPairList);
		static CubStorageBuffer AllocateCubHistogramBuffer(VolumeEdepPair edepPairList, int* histogramVals, double* logBins, int nbins);
		
		void* storage;
		size_t size;
};