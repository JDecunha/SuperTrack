//SuperTrack
#include "SuperTrackManager.hh"
#include "SimulationMethodFactory.hh"
#include "Track.cuh"
#include "ThreadAllocator.hh"
//ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TEntryList.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TMath.h"
#include "ROOT/TProcessExecutor.hxx"
#include "TCanvas.h"
#include "THStack.h"
#include "TPad.h"
//STD
#include <iostream>
#include <sstream>

void SuperTrackManager::AddSimulationMethod(const std::string& name, SimulationMethod* (*constructorPointer)(const INIReader&))
{
	//Add the name and simulation method pointer to the available factory constructors
	SimulationMethodFactory::GetInstance().AddSimulationMethod(name, constructorPointer);
}

void SuperTrackManager::Initialize(INIReader* reader)
{
	_inputFileReader = reader; //Store the INIReader pointer locally
	_output.SetDirectory(nullptr); //Have the ROOT histogram not be associated with the gROOT

	//Initialize the thread allocations with information from the reader
	InitializeThreadAllocations();

	//Set ROOT to be thread safe
	ROOT::EnableThreadSafety();

	//If thread allocation creation is succesful then we consider the initialization complete
	_bInitialized = true;
}

void SuperTrackManager::InitializeThreadAllocations()
{
	ThreadAllocator folderAllocator = ThreadAllocator(*_inputFileReader);
	folderAllocator.ReturnThreadAllocations(_threadAllocations);
}

void SuperTrackManager::Run()
{
	if (_bInitialized)
	{
		auto threadProcess = [=](ThreadAllocation threadInput)
		{
			Histogram histogram = Histogram(*_inputFileReader); //Construct thread local histogram
			SimulationMethod* method = SimulationMethodFactory::GetInstance().Construct(*_inputFileReader); //Construct thread local simulation method
			
			for (ThreadTask task : threadInput.GetTasks()) //Loop through each task for this thread
			{
				//Allocate memory for and load the track on GPU
				Track track; track.AllocateAndLoadTrack(task);
				
				//Allocate GPU only memory for the volume:edep paired list
				VolumeEdepPair edepsInTarget; edepsInTarget.Allocate(task.GetExitPoint()-task.GetEntryPoint());

				//Allocate required memory in method and histogram, for this track
				method->AllocateTrackProcess(track, task);
				histogram.AllocateTrackProcess(edepsInTarget);

				//Process the track oversample number of times
				for (int oversample = 0; oversample < task.GetNOversamples(); oversample++)
				{
					method->ProcessTrack(track, edepsInTarget); //Take the track, and return edepsInTarget
					histogram.SortReduceAndAddToHistogram(edepsInTarget); //From edepInTarget, add to the histogram
				}
				cudaDeviceSynchronize();
				//Free only memory allocated during track processing
				method->FreeTrackProcess();
				histogram.FreeTrackProcess(); 
				//TODO: have track and volumeedeppair destructors take care of this
				track.Free();
				edepsInTarget.Free();
			}

			//Clear the heap allocated method. Stack allocated class destructors are taken care of automatically
			delete method;
			
			return histogram.GetCPUHistogram();
		};

		// Create the pool of workers, based on the number of threads requested
		ROOT::TProcessExecutor workers(_threadAllocations.size());
		//Process the jobs and get a vector of the output
		std::vector<TH1D> process_output = workers.Map(threadProcess, _threadAllocations);
		//Reduce the output from each thread
		_output = Histogram::ReduceVector(process_output); //process output is left in an undefined state after this function

		//Export the final histogram
		EndOfRun(_output);
	}
	else {std::cout << "Can not start run without initializing manager." << std::endl;}
}

void SuperTrackManager::EndOfRun(TH1D& output)
{
	//Parse the information from the .ini
	std::string output_folder = _inputFileReader->Get("Output","Path","");
	std::string output_prefix = _inputFileReader->Get("Output","NamePrefix","");
	//Grab the random seed
	Long_t random_seed = _threadAllocations[0].GetRandomSeed();

	//Concatenate
	std::stringstream concatstream;
	concatstream << output_folder << output_prefix << random_seed << ".root";
	TString filename = concatstream.str();

	//Save
	TFile outfile(filename,"RECREATE");
	output.Write();

	//Output the number of tracks analyzed
	double nTracksAnalyzed = 0;
	int nOversamples = std::stoi(_inputFileReader->Get("Run","Oversamples",""));

	for (auto & allocation:_threadAllocations)
	{
		nTracksAnalyzed += allocation.GetTasks().size();
	}

	 TNamed tnNumberOfTracks = TNamed("Number of Tracks",std::to_string(nTracksAnalyzed));
	 TNamed tnNumberOfOversamples = TNamed("Number of Oversamples",std::to_string(nOversamples));
	 TNamed tnEffectiveNumberOfTracks = TNamed("Effective Number of Tracks",std::to_string(nTracksAnalyzed*nOversamples));

	 tnNumberOfTracks.Write();
	 tnNumberOfOversamples.Write();
	 tnEffectiveNumberOfTracks.Write();
}