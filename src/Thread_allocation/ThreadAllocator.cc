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
#include <algorithm>
#include <filesystem>
#include <iostream>


ThreadAllocator::ThreadAllocator(const std::string& folderPath, const int& numThreads, const int& nOversamples, const int& lowerFileLimit, const int& upperFileLimit, const Long_t& randomSeed)
{
	_folderPath = folderPath;
	_numThreads = numThreads;
	_lowerFileLimit = lowerFileLimit;
	_upperFileLimit = upperFileLimit;
	_randomSeed = randomSeed;
	_nOversamples = nOversamples;
}

ThreadAllocator::ThreadAllocator(const INIReader& reader)
{
	//Check that analysis over a folder is requested
	if (reader.Get("Run","Type","") == "Folder")
	{
		//Parse the file, return error if no file named
		_folderPath = reader.Get("Run","Path","");
		if (_folderPath == "") {std::cout << "Macro file error: ([Run], Folder) not defined." << std::endl; abort();}

		//Parse the other inputs. These are optional, no error if these are not present.
		_lowerFileLimit = reader.GetInteger("Run","FirstFile",-1);
		_upperFileLimit = reader.GetInteger("Run","LastFile",-1);
		_numThreads = reader.GetInteger("Run","Threads",1);
		_nOversamples = reader.GetInteger("Run","Oversamples",1);

		//Check if random seed argument present, if not, take the current time
		_randomSeed = reader.GetInteger("Run","RandomSeed",-1);
		if (_randomSeed == -1) { _randomSeed = time(NULL); }
	}
	else {std::cout << "Macro file error: ([Run], Type) does not = folder. Only analysis runs over folders are supported." << std::endl; abort();}
}

void ThreadAllocator::ReturnThreadAllocations(std::vector<ThreadAllocation>& threadAllocations)
{
	//Part 1: Save paths of .root files in the directory
	std::vector<std::string> filePaths;

	for(const auto &file : std::filesystem::directory_iterator(_folderPath))
		if(file.path().extension() == ".root")
			filePaths.push_back(file.path());
		
	//Sort the paths alphabetically
	std::sort(filePaths.begin(),filePaths.end());
	
	//Part 2: parse the lower and upper file limit arguments
	ParseFileLimits(filePaths.size());

	//Part 3: Mask a task for every track in the requested files
	std::vector<ThreadTask> tasks;

	for (int i = _lowerFileLimit; i <= _upperFileLimit; i++)
		//i-1, because _lowerFileLimit numbering starts at 1
		MakeTasks((TString)(filePaths[i-1]), tasks); //Get the tasks from the current file
	

	//Modify the number of threads allocated to if there are fewer tracks than threads
	if (tasks.size() < _numThreads)
	{
		std::cout << "There are fewer tracks to analyze than the number of requested threads." << std::endl;
		std::cout << "Lowering number of active threads to: " << tasks.size() << " to accomodate." << std::endl;
		_numThreads = tasks.size();
	}

	//Part 4: Create ThreadAllocations (which are a collection of tasks)
	//for every thread. Depending on # of tracks and threads
	int numTracksPerThread = tasks.size()/_numThreads; //this is how many tracks we allocate per thread
	int numTracksPerThreadRemainder = tasks.size() % _numThreads;  //and for remainder# of threads we allocate an additional track
	int currentTrack = 0; //this is used in the loop below

	// *** Loop logic explanation: ***
	//Imagine 15 cores analyzing 20 tracks
	//the numTracksPerThread is 1, remainder 5
	//First 5 cores analyze numTracksPerThread+1 tracks, 2 tracks per thread
	//Last 10 cores just analyze numTracksPerThread tracks, 1 track per thread
	//Works out to 20 tracks

	for (int i = 0; i < _numThreads; i++) //loop for each thread
	{
		threadAllocations.push_back(ThreadAllocation(i,_randomSeed,_nOversamples)); //create the allocation for a thread

		int numTasksToCreate = numTracksPerThread;

		if (i < numTracksPerThreadRemainder) //add an extra task/track for these threads
			numTasksToCreate += 1;

		for (int j = 0; j < numTasksToCreate; j++) //create tasks and push them into the ThreadAllocation 
		{
			threadAllocations[i].AddTask(tasks[currentTrack]);
			currentTrack += 1;
		}
	}
}

void ThreadAllocator::ParseFileLimits(const int& folderSize)
{
	if (_lowerFileLimit == -1 && _upperFileLimit == -1) //default case, analyze whole folder
	{ 
		_lowerFileLimit = 1; //Numbering starts from 1
		_upperFileLimit = folderSize; 
	}
	else if (_upperFileLimit < _lowerFileLimit) //This case is invalid
	{
		std::cout << "Lower file limit greater than upper file limit. Aborting." << std::endl;
		abort();
	}
	else if (_upperFileLimit > folderSize)
	{
		std::cout << "Upper file limit greater than number of files. Aborting." << std::endl;
		abort();
	}
	else if (_lowerFileLimit < 1)
	{
		std::cout << "Lower file limit less than 1. Aborting." << std::endl;
		abort();
	}
}

void ThreadAllocator::MakeTasks(const TString& file, std::vector<ThreadTask>& inputTasks)
{
	//open the file and retrieve the number of tracks
	TFile f = TFile(file);
	TTree *trackIndex;
	f.GetObject("Track index",trackIndex);
	int nTracksInFile = trackIndex->GetEntries();

	//Establish a pointer to and readout the end point of the tracks
	long start_entry_val = 0;
	TTreeReader trackIndexReader("Track index", &f);
	TTreeReaderValue<long long> end_entry_val(trackIndexReader, "index");

	//Make the output vector
	std::vector<ThreadTask> output;

	for (Int_t i = 0; i < nTracksInFile; i++)
	{
		trackIndexReader.Next();
		inputTasks.push_back(ThreadTask(file,start_entry_val,*end_entry_val));
		start_entry_val = *end_entry_val;
	}
}