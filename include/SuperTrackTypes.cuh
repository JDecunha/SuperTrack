#pragma once

//ROOT
#include "TMath.h"
//inih
#include "INIReader.h"
//CUB
#include <cub/cub.cuh>

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

	SphericalGeometry(INIReader inputReader)
	{
		//Setting values
		scoringRegionHalfLength = inputReader.GetReal("VoxelConstrainedSphere","ScoringRegionHalfLength",0);
		sphereDiameter = inputReader.GetReal("VoxelConstrainedSphere","ScoringSphereDiameter",0);
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

struct CUBAddOperator
{
    template <typename T>
    CUB_RUNTIME_FUNCTION __forceinline__
    T operator()(const T &a, const T &b) const {
        return a+b;
    }
};