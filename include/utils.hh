#pragma once

//The compiler needs templated types to include their whole definition
//I can't just forward declare the classes I believe
#include "TH1F.h"
#include "Math/SMatrix.h"
typedef ROOT::Math::SMatrix<Double_t,3> SMatrix33;

//A header file for defining various utility functions
namespace CPUHistogramUtils
{
	//Histogram utility functions
	void BinLogX(TH1* h);
	void BinLogXMultithread(std::shared_ptr<TH1> h);
	void PMF_to_PDF(TH1* h);
	void Prepare_for_Semilog(TH1* h);
};

namespace utils
{
	void LogSpace(float bottom_order_mag, float top_order_mag, int nbins, double* bins);
};

//Other utility functions
void uniform_random_rotation_matrix_optimized(float x0,float x1,float x2, SMatrix33* matrix);