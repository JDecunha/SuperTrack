#include "SuperTrackManager.hh"
#include "SimulationMethodFactory.hh"
#include "Track.cuh"
#include "ThreadAllocator.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TEntryList.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TMath.h"
#include "ROOT/TProcessExecutor.hxx"

#include <iostream>

void SuperTrackManager::AddSimulationMethod(const std::string& name, SimulationMethod* (*constructorPointer)(const INIReader&))
{
	//Add the name and simulation method pointer to the available factory constructors
	SimulationMethodFactory::GetInstance().AddSimulationMethod(name, constructorPointer);
}

void SuperTrackManager::Initialize(INIReader* reader)
{
	_inputFileReader = reader; //Store the INIReader poitner locally

	//Initialize the thread allocations with information from the reader
	InitializeThreadAllocations();

	//If thread allocation creation is succesful then we consider the initialization complete
	_bInitialized = true;
}

void SuperTrackManager::InitializeThreadAllocations()
{
	if (_inputFileReader->Get("Run","Type","") == "Folder")
	{
		//Parse the file
		std::string folderPath = _inputFileReader->Get("Run","Path","");
		int firstFile = _inputFileReader->GetInteger("Run","FirstFile",-1);
		int lastFile = _inputFileReader->GetInteger("Run","LastFile",-1);
		int numThreads = _inputFileReader->GetInteger("Run","Threads",1);
		int numOversamples = _inputFileReader->GetInteger("Run","Oversamples",1);

		//Create the thread allocator
		ThreadAllocator folderAllocator = ThreadAllocator(folderPath,numThreads,numOversamples,firstFile,lastFile);

		//Store the allocations in this class
		folderAllocator.ReturnThreadAllocations(_threadAllocations);
	}
	else {std::cout << "Only analysis runs over folders are supported." << std::endl; abort();}
}


void SuperTrackManager::Run()
{
	if (_bInitialized)
	{
		auto threadProcess = [=](ThreadAllocation threadInput)
		{
			//Retrieve the vector of tasks
			std::vector<ThreadTask> tasks = threadInput.GetTasks();

			//Construct thread local histograms
			Histogram histogram = Histogram(*_inputFileReader);

			//Construct thread local simulation method
			SimulationMethodFactory& methodFactory = SimulationMethodFactory::GetInstance();
			SimulationMethod* method = methodFactory.Construct(*_inputFileReader);

			//For each task in the ThreadAllocation
			for (ThreadTask task : tasks)
			{
				//Load the track
				Track track;
				track.AllocateAndLoadTrack(task);
				
				//Allocate GPU only memory for the volume:edep paired list
				VolumeEdepPair edepsInTarget;
				edepsInTarget.Allocate(task.GetExitPoint()-task.GetEntryPoint()+1);

				//Allocate required memory for simulation method
				method->AllocateTrackProcess(track,task);

				//Allocate required buffer memory for the histogram
				histogram.Allocate(edepsInTarget);

				for (int oversample = 0; oversample < task.GetNOversamples(); oversample++)
				{
					method->ProcessTrack(track,edepsInTarget);
					histogram.SortReduceAndAddToHistogram(edepsInTarget);
				}

				int number_of_values_in_histogram = 0;
				cudaDeviceSynchronize();
				//Read out histogram
				for (int i = 0; i < histogram._nbins; i++)
				{
					number_of_values_in_histogram += histogram._histogramValsAccumulated[i];
					std::cout << "Bin: " << histogram._binEdges[i] << " Counts: " << histogram._histogramValsAccumulated[i] << std::endl;
				}
				std::cout << number_of_values_in_histogram << std::endl;

				method->Free();
				histogram.Free(); 
			}

			TH1F lineal_histogram = TH1F("Lineal energy histogram", "y*f(y)", 200, -2,1);
			return lineal_histogram;
		};

		// Create the pool of workers
		ROOT::TProcessExecutor workers(_threadAllocations.size());
		//Process the jobs and get a vector of the output
		std::vector<TH1F> process_output = workers.Map(threadProcess, _threadAllocations);


	}
	else {std::cout << "Can not start run without initializing manager." << std::endl;}
}
