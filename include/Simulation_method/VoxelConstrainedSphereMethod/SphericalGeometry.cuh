#pragma once

//ROOT
#include "TMath.h"
//inih
#include "INIReader.h"

struct SphericalGeometry
{
	SphericalGeometry(double scoring_region_half_length, double sphere_diameter);
	SphericalGeometry(INIReader inputReader);

	//The geometry
	double scoringRegionLength;
	double scoringRegionHalfLength;
	double greatestSphereOffset;
	int numSpheresLinear;
	
	//The targets
	double sphereDiameter;
	double sphereRadius;
	
};