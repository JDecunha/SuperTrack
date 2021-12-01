#pragma once

//ROOT
#include "TMath.h"
//inih
#include "INIReader.h"

struct SphericalGeometry
{
	SphericalGeometry(double scoring_region_half_length, double sphere_diameter);
	SphericalGeometry(INIReader inputReader);

	double scoringRegionHalfLength;
	double greatestSphereOffset;
	double sphereDiameter;
	double sphereRadius;
	int numSpheresLinear;
};