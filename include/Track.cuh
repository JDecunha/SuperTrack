#pragma once

#include <tuple>
#include "TROOT.h"

class Track
{ 
	public:

		Track(int nVals);
		~Track();
		static void LoadTrack(const std::tuple<Int_t,Int_t,Int_t,TString> &input, Track *deviceTrack);

		double* x;
		double* y;
		double* z;
		double* edep; 
};

