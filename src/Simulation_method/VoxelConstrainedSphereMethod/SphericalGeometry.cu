#include "SphericalGeometry.cuh"

//This is just a helper class for VoxelConstrainedSpereMethod
//Spherical geometry structs get sent to the GPU kernels

//define the two constructors
SphericalGeometry::SphericalGeometry(double scoring_region_half_length, double sphere_diameter)
{
	//Setting values
	scoringRegionHalfLength = scoring_region_half_length;
	sphereDiameter = sphere_diameter;
	sphereRadius = sphereDiameter/2;

	//The number of spheres in a line set with TMath::Floor, so that the box can fit all targets without cutting any off
	numSpheresLinear = TMath::Floor(((scoringRegionHalfLength*2)/sphereDiameter)); 
	
	//We scale the size of the box to match the number of spheres that fit within it
	scoringRegionLength = numSpheresLinear*sphereDiameter;
	scoringRegionHalfLength = scoringRegionLength/2;

	//gSO is always just one radius off the edge
	greatestSphereOffset = -scoringRegionHalfLength+sphereRadius;

	//Check that radius is not larger than half length
	if(sphereRadius > scoringRegionHalfLength)
	{
		scoringRegionHalfLength = sphereRadius;
	}
}

SphericalGeometry::SphericalGeometry(INIReader inputReader)
{
	//Setting values from the .ini file
	scoringRegionHalfLength = inputReader.GetReal("VoxelConstrainedSphere","ScoringRegionHalfLength",0);
	sphereDiameter = inputReader.GetReal("VoxelConstrainedSphere","ScoringSphereDiameter",0);
	sphereRadius = sphereDiameter/2;

	//The number of spheres in a line set with TMath::Floor, so that the box can fit all targets without cutting any off
	numSpheresLinear = TMath::Floor(((scoringRegionHalfLength*2)/sphereDiameter)); 
	
	//We scale the size of the box to match the number of spheres that fit within it
	scoringRegionLength = numSpheresLinear*sphereDiameter;
	scoringRegionHalfLength = scoringRegionLength/2;

	//gSO is always just one radius off the edge
	greatestSphereOffset = -scoringRegionHalfLength+sphereRadius;

	//Check that radius is not larger than half length
	if(sphereRadius > scoringRegionHalfLength)
	{
		scoringRegionHalfLength = sphereRadius;
	}
}
