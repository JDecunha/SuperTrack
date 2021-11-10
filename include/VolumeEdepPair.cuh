#pragma once

#include "stdint.h"

class VolumeEdepPair
{
	public:
		VolumeEdepPair();

		void Allocate(uint64_t numInputElements);
		void Free();

		uint64_t* volume;
		double* edep;
		int* numElements; //this is type pointer but it should only point to a single value
};
