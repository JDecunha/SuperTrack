#pragma once

#include "SimulationMethod.hh"
#include "TROOT.h"
#include "SuperTrackTypes.cuh"
#include "Track.cuh"


class VoxelConstrainedSphereMethod : public SimulationMethod
{
	public:
		VoxelConstrainedSphereMethod(const INIReader& macroReader);

		void ParseInput() override;
		void AllocateTrackProcess(Track track, ThreadTask task) override; //This will handle all the memory allocations
		VolumeEdepPair ProcessTrack(Track track) override;
		void FreeTrackProcess() override;
		void Free() override;

	private:
		Track _randomlyShiftedTrack;
		//SphericalGeometry _sphericalGeometry;
		int* _numInBox;
		int* _inSphereTrackId;
		float* _randomVals;

};

namespace VoxelConstrainedSphereMethodKernel
{
	__global__ void ScoreTrackInSphere(SphericalGeometry geometry, Track inputTrack, int *numElements, int *trackIdInSphere, VolumeEdepPair outputPair);
	__global__ void FilterTrackInSphere(SphericalGeometry geometry, Track inputTrack, int *numElements, int *numElementsCompacted, int *trackIdInSphere);
	__global__ void FilterInScoringBox(SphericalGeometry geometry, float* randomVals, Track inputTrack, Track outputTrack, int numElements, int *numElementsCompacted, int oversampleIterationNumber);
};