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
//INIH
#include "INIReader.h"

//The job of ThreadAllocator:
//is to take the input arguments regarding which folder/files to analyze
//And split them up so they can be run on many cores

class ThreadAllocator
{
	public:
		ThreadAllocator(const std::string& folderPath, const int& numThreads, const int& nOversamples = 1, const int& lowerFileLimit = -1, const int& upperFileLimit = -1, const Long_t& randomSeed = time(NULL));

		ThreadAllocator(const INIReader& reader);

		void ReturnThreadAllocations(std::vector<ThreadAllocation>& threadAllocations);

	private:
		void MakeTasks(const TString& file, std::vector<ThreadTask>& inputTasks);
		void ParseFileLimits(const int& folderSize); //Take the lower and upper file limits, and error check

		std::string _folderPath;
		int _numThreads;
		int _nOversamples;
		int _lowerFileLimit;
		int _upperFileLimit;
		long _randomSeed;
};