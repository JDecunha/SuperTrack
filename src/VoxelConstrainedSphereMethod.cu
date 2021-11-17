#include "VoxelConstrainedSphereMethod.hh"
#include "SuperTrackTypes.cuh"
#include "VolumeEdepPair.cuh"


VoxelConstrainedSphereMethod::VoxelConstrainedSphereMethod(const INIReader& macroReader) : SimulationMethod(macroReader)
{
	ParseInput();
}

//ParseInput takes the INIReader and parses the file to initialize the class
void VoxelConstrainedSphereMethod::ParseInput()
{
	double scoringRegionHalfLength = _macroReader.GetReal("VoxelConstrainedSphere","ScoringRegionHalfLength",0);
	double scoringSphereDiameter = _macroReader.GetReal("VoxelConstrainedSphere","ScoringSphereDiameter",0);

	std::cout << scoringRegionHalfLength << std::endl;
	std::cout << scoringSphereDiameter << std::endl;
}

void VoxelConstrainedSphereMethod::AllocateTrackProcess(Track track, ThreadTask task) 
{ 

}

VolumeEdepPair VoxelConstrainedSphereMethod::ProcessTrack(Track track)
{ 
	VolumeEdepPair huh;
	return huh;
}

void VoxelConstrainedSphereMethod::FreeTrackProcess()
{ 

}

void VoxelConstrainedSphereMethod::Free()
{ 

}

