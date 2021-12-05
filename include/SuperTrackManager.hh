#pragma once

#include "Histogram.cuh"
#include "SimulationMethod.hh"
#include "ThreadAllocation.hh"
#include "INIReader.h"
#include <vector>

class SuperTrackManager
{	
	public:
		//Proper method to instantiate singleton
		static SuperTrackManager& GetInstance()
		{
			static SuperTrackManager singleton;

			return singleton;
		}

		~SuperTrackManager() { delete _inputFileReader; }

		//Delete the copy and = constructors
		SuperTrackManager(SuperTrackManager const&) = delete;
		void operator=(SuperTrackManager const&) = delete;

		//Will access the method factory and add a new method
		void AddSimulationMethod(const std::string& name, SimulationMethod* (*constructorPointer)(const INIReader&));

		//Will take and store the INIReader pointer, and populate thread allocations
		void Initialize(INIReader* reader);

		//Passes the INIreader and thread allocations to each thread
		//Creates thread local histograms and simulation method
		void Run();

	private:
		//private singleton constructor
		SuperTrackManager() {_bInitialized = false; } 

		void InitializeThreadAllocations();
		void EndOfRun(TH1D& output);
		
		//Internal values
		INIReader* _inputFileReader;
		std::vector<ThreadAllocation> _threadAllocations;
		bool _bInitialized;
};