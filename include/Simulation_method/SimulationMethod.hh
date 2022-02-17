#pragma once

class Track;
class VolumeEdepPair;
#include "INIReader.h"
#include "ThreadTask.hh"
//ROOT
#include "TROOT.h"

//So a concrete version of a SimulationMethod is a combined geometry AND random
//shift generating class. Those two things together specify
//a simulation method.

//What is defined here is just an interface.

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