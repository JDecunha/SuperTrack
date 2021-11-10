#pragma once

#include "TMath.h"
#include <cub/cub.cuh>

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

//Miscellaneous

struct CUBAddOperator
{
    template <typename T>
    CUB_RUNTIME_FUNCTION __forceinline__
    T operator()(const T &a, const T &b) const {
        return a+b;
    }
};