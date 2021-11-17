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
//I could call it SuperimposeMethod too.

class SimulationMethod
{
	public:
		SimulationMethod(const INIReader& macroReader);
		
		virtual void ParseInput() = 0;
		virtual void AllocateTrackProcess(Track track, ThreadTask task) = 0; //This will handle all the memory allocations
		virtual VolumeEdepPair ProcessTrack(Track track) = 0;
		virtual void FreeTrackProcess() = 0;
		virtual void Free() = 0;

	protected:
		INIReader _macroReader;
};

namespace SimulationMethodKernel
{
	__global__ void ZeroInt(int* toZero);
	void GenerateRandomXYShift(const std::tuple<Int_t,Int_t,Int_t,TString> &input, float **randomVals, const int &nSamples, const long &random_seed);

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