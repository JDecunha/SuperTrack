#pragma once

#include <tuple>
#include "TROOT.h"
class ThreadTask;

class Track
{ 
	public:

		Track();
		void Free();
		void AllocateAndLoadTrack(const std::tuple<Int_t,Int_t,Int_t,TString>& input);
		void AllocateAndLoadTrack(const ThreadTask& task);
		void AllocateEmptyTrack(int nVals);


		double* x;
		double* y;
		double* z;
		double* edep; 
};

