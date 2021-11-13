//The job of ThreadFileAllocator:
//is to take the input arguments regarding which folder/files to analyze
//And split them up so they can be run on many cores

//SuperTrack
#include "ThreadAllocator.hh"
//ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TMath.h"
//STD library
#include <filesystem>
#include <iostream>


ThreadAllocator::ThreadAllocator(const std::string& folderPath, const int& numThreads, const int& lowerFileLimit, const int& upperFileLimit, const Long_t& randomSeed)
{
	_folderPath = folderPath;
	_numThreads = numThreads;
	_lowerFileLimit = lowerFileLimit;
	_upperFileLimit = upperFileLimit;
	_randomSeed = randomSeed;
}

std::vector<ThreadAllocation> ThreadAllocator::ReturnThreadAllocations()
{
	//Part 1: determine number of .root files and save their paths
	std::vector<std::string> filePaths;
	int numFiles = 0;

	for(const auto &file : std::filesystem::directory_iterator(_folderPath))
	{ 
		if(file.path().extension() == ".root")
		{
			filePaths.push_back(file.path());
		}
	}
	
	//Part 2: parse the lower and upper file limit arguments
	if (_lowerFileLimit == 0 && _upperFileLimit == 0) //default case, analyze whole folder
	{ 
		_lowerFileLimit = 1; //Numbering starts from 1 for this
		_upperFileLimit = filePaths.size(); 
	}
	else if (_upperFileLimit != 0 && _upperFileLimit < _lowerFileLimit) //This case is invalid
	{
		std::cout << "Lower file limit greater than upper file limit. Aborting." << std::endl;
		abort();
	}

	//Part 3: for the requested files (within the limits)
	//make a task for every Track in the file
	std::vector<ThreadTask> tasks;
	int numTracks;
	int numTracksTotal = 0;

	//Make a task out of every track in the requested files
	for (int i = _lowerFileLimit; i <= _upperFileLimit; i++)
	{
		numTracks = MakeTasks((TString)(filePaths[i-1]), tasks); //Get the tasks from the current file
		numTracksTotal += numTracks;
	}

	//Abort in this case
	if (numTracksTotal < _numThreads)
	{
		std::cout << "There are fewer tracks to analyze than the number of requested threads." << std::endl;
		std::cout << "Lowering number of active threads to: " << numTracksTotal << " to accomodate." << std::endl;
		_numThreads = numTracksTotal;
	}

	//Part 4: Create out ThreadAllocations (which are a collection of tasks)
	//for every thread. Depending on # of tracks and threads

	std::vector<ThreadAllocation> threadAllocations; //Create the output thread allocations
	int numTracksPerThread = 0; //this is how many tracks we allocate per thread
	int numTracksPerThreadRemainder = 0; //and for remainder# of threads we allocate ONE more track
	int currentTrack = 0;

	// *** Loop logic explanation: ***
	//Imagine 10 cores analyzing 15 files
	//the numTracksPerThread is 1, remainder 5
	//First 5 cores analyze numTracksPerThread+1 tracks, 2 tracks per thread
	//Last 5 cores just analyze numTracksPerThread tracks, 1 track per thread
	//Works out to 15 tracks


	//Get the quotient and remainder
	numTracksPerThread = numTracksTotal/_numThreads;
	numTracksPerThreadRemainder = numTracksTotal % _numThreads; //Remainder exists if numTracks not divisible by numCores

	for (int i = 0; i < _numThreads; i++) //loop over each requested thread
	{
		int j = currentTrack;
		ThreadAllocation allocation = ThreadAllocation(_randomSeed,i); //create the allocation for a thread

		if (i < numTracksPerThreadRemainder) //add an extra track for these threads
		{
			while (currentTrack < (j+numTracksPerThread+1)) //Loop over and add tasks to the allocation
			{
				allocation.AddTask(tasks[currentTrack]);
				currentTrack += 1;
			}
		}
		else //don't add an extra track for these
		{
			while (currentTrack < (j+numTracksPerThread)) //Loop over and add tasks to the allocation
			{
				allocation.AddTask(tasks[currentTrack]);
				currentTrack += 1;
			}
		}

		//add the allocation for this thread to the global list
		threadAllocations.push_back(allocation); 
	}

	return threadAllocations;
}

int ThreadAllocator::GetNumberOfTracks(const TString& file)
{
	//open the file and retrieve the number of tracks
	TFile f = TFile(file);
	TTree *trackIndex;
	f.GetObject("Track index",trackIndex);
	int nTracksInFile = trackIndex->GetEntries();

	return nTracksInFile;
}

int ThreadAllocator::MakeTasks(const TString& file, std::vector<ThreadTask>& inputTasks)
{
	//open the file and retrieve the number of tracks
	TFile f = TFile(file);
	TTree *trackIndex;
	f.GetObject("Track index",trackIndex);
	int nTracksInFile = trackIndex->GetEntries();

	//Establish a pointer to and readout the tracks
	long start_entry_val = 0;
	TTreeReader trackIndexReader("Track index", &f);
	TTreeReaderValue<long long> end_entry_val(trackIndexReader, "index");

	//Make the output vector
	std::vector<ThreadTask> output;

	for (Int_t i = 0; i < nTracksInFile; i++)
	{
		trackIndexReader.Next();
		output.push_back(ThreadTask(file,start_entry_val,*end_entry_val));
		start_entry_val = *end_entry_val;
	}

	inputTasks.insert(inputTasks.end(),output.begin(),output.end()); //Append them to the global vector of tasks

	return nTracksInFile;
}