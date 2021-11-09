#pragma once

#include "TMath.h"

//
//Structs and typedefs
//

//Geometry

struct SphericalGeometry
{
	SphericalGeometry(double scoring_region_half_length, double sphere_diameter)
	{
		//Setting values
		scoringRegionHalfLength = scoring_region_half_length;
		sphereDiameter = sphere_diameter;
		sphereRadius = sphereDiameter/2;

		//Calculating values
		numSpheresLinear = TMath::Ceil(((scoringRegionHalfLength*2)/sphereDiameter)); 
		greatestSphereOffset = -(((float(numSpheresLinear))/2)-0.5)*sphereDiameter;
	}

	double scoringRegionHalfLength;
	double greatestSphereOffset;
	double sphereDiameter;
	double sphereRadius;
	int numSpheresLinear;
};

//Track-related

/*struct Track
{ 
	double x;
	double y;
	double z;
	double edep; 
};*/

struct Track
{ 
	double* x;
	double* y;
	double* z;
	double* edep; 
};

struct VolumeEdepPair
{
	uint64_t* volume;
	double* edep;
	int* numElements; //this is type pointer but it should only point to a single value
};

//Miscellaneous

struct CubStorageBuffer
{
	CubStorageBuffer()
	{
		storage = NULL;
		size = 0;
	}
	void* storage;
	size_t size;
};