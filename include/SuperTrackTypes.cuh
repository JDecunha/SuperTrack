#pragma once

#include "TMath.h"

//
//Structs and typedefs
//

struct Track
{ 
	double x;
	double y;
	double z;
	double edep; 
};
typedef struct Track Track; //this is a typdef which maps from struct Track --> Track. Just saves us from writing struct when we refer to the struct later.

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
typedef struct SphericalGeometry SphericalGeometry;

struct VolumeEdepPair
{
	uint64_t* volume;
	double* edep;
	int* numElements; //this is type pointer but it should only point to a single value
};
typedef struct VolumeEdepPair VolumeEdepPair;