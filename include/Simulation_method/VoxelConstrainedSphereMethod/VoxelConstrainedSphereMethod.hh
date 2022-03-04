#pragma once

#include "SimulationMethod.hh"
#include "TROOT.h"
#include "SphericalGeometry.cuh"
#include "Track.cuh"


class VoxelConstrainedSphereMethod : public SimulationMethod
{
	public:
		VoxelConstrainedSphereMethod(const INIReader& macroReader);

		void ParseInput() override;
		void AllocateTrackProcess(Track track, ThreadTask task) override; //This will handle all the memory allocations
		void ProcessTrack(Track track, VolumeEdepPair& edepsInTarget) override;
		void FreeTrackProcess() override;
		void Free() override;

		static SimulationMethod* Construct(const INIReader& macroReader);

	private:

		//Number of threads and blocks suggested by user
		int _suggestedCudaThreads;
		int _suggestedCudaBlocks;
		bool _randomShifts;

		//Takes the information from the .ini reader and make a geometry
		SphericalGeometry _sphericalGeometry;

		//These values are allocated before a track is processed and used when superimposing the Track
		int _nSteps;
		int _oversampleIterationNumber;
		float* _randomVals;
		Track _randomlyShiftedTrack;
		int* _numInVoxel;
		int* _inSphereTrackId;
};

namespace VoxelConstrainedSphereMethodKernel
{
	__global__ void ScoreTrackInSphere(SphericalGeometry geometry, Track inputTrack, int *numElements, int *trackIdInSphere, VolumeEdepPair outputPair);
	__global__ void FilterTrackInSphere(SphericalGeometry geometry, Track inputTrack, int *numElements, int *numElementsCompacted, int *trackIdInSphere);
	__global__ void FilterInScoringBox(SphericalGeometry geometry, float* randomVals, Track inputTrack, Track outputTrack, int numElements, int *numElementsCompacted, int oversampleIterationNumber);
};