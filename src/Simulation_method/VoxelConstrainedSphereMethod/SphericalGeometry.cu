#include "SphericalGeometry.cuh"

//define the two constructors
//This is just a helper class for VoxelConstrainedSpereMethod
//Spherical geometry structs get sent to the GPU kernels

SphericalGeometry::SphericalGeometry(double scoring_region_half_length, double sphere_diameter)
{
	//Setting values
	scoringRegionHalfLength = scoring_region_half_length;
	sphereDiameter = sphere_diameter;
	sphereRadius = sphereDiameter/2;

	//Calculating values
	numSpheresLinear = TMath::Ceil(((scoringRegionHalfLength*2)/sphereDiameter)); 
	greatestSphereOffset = -(((float(numSpheresLinear))/2)-0.5)*sphereDiameter;
}

SphericalGeometry::SphericalGeometry(INIReader inputReader)
{
	//Setting values
	scoringRegionHalfLength = inputReader.GetReal("VoxelConstrainedSphere","ScoringRegionHalfLength",0);
	sphereDiameter = inputReader.GetReal("VoxelConstrainedSphere","ScoringSphereDiameter",0);
	sphereRadius = sphereDiameter/2;

	//Calculating values
	numSpheresLinear = TMath::Ceil(((scoringRegionHalfLength*2)/sphereDiameter)); 
	greatestSphereOffset = -(((float(numSpheresLinear))/2)-0.5)*sphereDiameter;
}