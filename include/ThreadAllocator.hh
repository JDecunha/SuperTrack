#pragma once

//SuperTrack
#include "ThreadAllocation.hh"
#include "ThreadTask.hh"
//STD
#include <tuple>
#include <vector>
#include "time.h"
//ROOT
#include "TROOT.h"

//The job of ThreadAllocator:
//is to take the input arguments regarding which folder/files to analyze
//And split them up so they can be run on many cores

class ThreadAllocator
{
	public:
		ThreadAllocator(const std::string &folderPath, const int &numThreads, const int &lowerFileLimit = 0, const int &upperFileLimit = 0, const Long_t& randomSeed = time(NULL));
		std::vector<ThreadAllocation> ReturnThreadAllocations();

	private:
		int GetNumberOfTracks(const TString& file);
		int MakeTasks(const TString& file, std::vector<ThreadTask>& inputTasks);

		std::string _folderPath;
		int _numThreads;
		int _lowerFileLimit;
		int _upperFileLimit;
		long _randomSeed;
};