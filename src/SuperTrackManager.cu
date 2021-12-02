#include "SuperTrackManager.hh"
#include "Track.cuh"

#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TEntryList.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TMath.h"
#include "ROOT/TProcessExecutor.hxx"

#include <iostream>

SuperTrackManager::SuperTrackManager() 
{
	_histogram = NULL;
	_method = NULL;
	_bInitialized = false;
}

void SuperTrackManager::AddThreadAllocations(std::vector<ThreadAllocation>& allocations)
{
	_threadAllocations = std::move(allocations);
}

void SuperTrackManager::AddHistogram(Histogram* histogram)
{
	_histogram = histogram;
}

void SuperTrackManager::AddSimulationMethod(SimulationMethod* method)
{
	_method = method;
}

void SuperTrackManager::Initialize()
{
	if (_histogram == NULL) { std::cout << "Manager can not initialize: does not possess valid histogram pointer." << std::endl; return; }
	if (_method == NULL) { std::cout << "Manager can not initialize: does not possess valid simulation method pointer." << std::endl; return; }
	if (_threadAllocations.empty()) { std::cout << "Manager can not initialize: does not have valid thread allocations." << std::endl; return; }

	_bInitialized = true;
} 

void SuperTrackManager::Run()
{
	if (_bInitialized)
	{
		auto threadProcess = [=](ThreadAllocation threadInput)
		{
			std::vector<ThreadTask> tasks = threadInput.GetTasks();

			//For each task in the ThreadAllocation
			for (ThreadTask task : tasks)
			{
				//Load the track
				Track deviceTrack;
				deviceTrack.AllocateAndLoadTrack(task);
				/*
				//Allocate GPU only memory for the volume:edep paired list
				VolumeEdepPair edepsInTarget;
				edepsInTarget.Allocate(task.GetExitPoint()-task.GetEntryPoint()+1);

				//Allocate required memory for simulation method
				_method->AllocateTrackProcess(deviceTrack,task);

				//Allocate required buffer memory for the histogram
				_histogram.Allocate(edepsInTarget);

				for (int oversampleIteration = 0; oversampleIteration < task.GetNOversamples(); oversampleIteration++)
				{
					_method.ProcessTrack(track,edepsInTarget);
					_histogram.SortReduceAndAddToHistogram(edepsInTarget);
				}
				cudaDeviceSynchronize();

				int number_of_values_in_histogram = 0;
				cudaDeviceSynchronize();
				//Read out histogram
				for (int i = 0; i < _histogram._nbins; i++)
				{
					number_of_values_in_histogram += _histogram._histogramValsAccumulated[i];
					std::cout << "Bin: " << _histogram._binEdges[i] << " Counts: " << _histogram._histogramValsAccumulated[i] << std::endl;
				}

				_method.Free();
				_histogram.Free(); */
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

void SuperTrackManager::TestSuperimpose(std::vector<ThreadAllocation>& allocations, Histogram* histogram, SimulationMethod* method)
{

}