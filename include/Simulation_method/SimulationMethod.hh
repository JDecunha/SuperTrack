#pragma once

class Track;
class VolumeEdepPair;
#include "INIReader.h"
#include "ThreadTask.hh"
//ROOT
#include "TROOT.h"

//So the simulationMethod is a combined geometry AND random
//shift generating class. Those two things together specify
//a simulation method.

class SimulationMethod
{
	public:
		SimulationMethod(const INIReader& macroReader);

		//Virtual functions
		virtual void ParseInput() = 0;
		virtual void AllocateTrackProcess(Track track, ThreadTask task) = 0; //This will handle all the memory allocations
		virtual void ProcessTrack(Track track, VolumeEdepPair& edepsInTarget) = 0;
		virtual void FreeTrackProcess() = 0;
		virtual void Free() = 0;

	protected:
		//Static helper functions for derived SimulationMethods
		static void GenerateRandomXYShift(const ThreadTask &task, float** randomVals);

		INIReader _macroReader;
};

namespace SimulationMethodKernel
{
	__global__ void ZeroInt(int* toZero);

};

//So the interface will be as follows:
//The SuperTrackRunManager will parse the input and look for the
//simulationMethod specifier
//Which will then prompt it to generate a simulationMethod of the
//specified type.
//
//The simulation method will receieve the config file as an input parameter
//allowing it to fill its required fields from the .ini
//If any field is missing from the .ini it will abort
//
//I have to think about how the seeding of random numbers will work
//as well here though
//Include another public function that gets called
//before ProcessTrack every time?


//What are my classes?
//Histogram class
//simulationMethod class
//Track class