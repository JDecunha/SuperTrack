#include "utils.hh"
#include "TMath.h"

void BinLogX(TH1F* h)
{
	///Very useful function from the ROOT forums
	///If you enter your axes in Log10 notation (i.e. -3,0) in the initial TH1 constructor
	///This will reformat your bins logarithmically for you

   TAxis *axis = h->GetXaxis();
   int bins = axis->GetNbins();

   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];

   for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);

   }
   axis->Set(bins, new_bins);
   delete[] new_bins;
} 

void BinLogXMultithread(std::shared_ptr<TH1F> h)
{
	///Very useful function from the ROOT forums
	///If you enter your axes in Log10 notation (i.e. -3,0) in the initial TH1 constructor
	///This will reformat your bins logarithmically for you

   TAxis *axis = h->GetXaxis();
   int bins = axis->GetNbins();

   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];

   for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);

   }
   axis->Set(bins, new_bins);
   delete[] new_bins;
} 

void PMF_to_PDF(TH1* h)
{
	//This converts a probability mass function to a probability density function
	int length = h->GetNbinsX();
	float_t normalization = 0;
	for (int i = 0; i < length; i++) //bins per bin width
	{
		auto value = h->GetBinContent(i);
		auto low_edge = h->GetBinLowEdge(i);
		auto high_edge = h->GetBinLowEdge(i+1);
		auto per_width = value/(high_edge-low_edge);
		h->SetBinContent(i,per_width);

		normalization += value;
	}
	for (int i = 0; i < length; i++) //normalize
	{
		auto value = h->GetBinContent(i);
		auto normalized = value/normalization;
		h->SetBinContent(i,normalized);
	}
}

void Prepare_for_Semilog(TH1* h)
{
	//This takes a PDF and normalizes it by an extra factor of y (i.e. to give y*f(y))to preserve the graphical properties
	//of a PDF on a semilog axis
	int length = h->GetNbinsX();
	for (int i = 0; i < length; i++) //bins per bin width
	{
		auto value = h->GetBinContent(i);
		auto low_edge = h->GetBinLowEdge(i);
		auto high_edge = h->GetBinLowEdge(i+1);
		auto mid_value = (high_edge+low_edge)/2;
		h->SetBinContent(i,value*mid_value);
	}
}

void uniform_random_rotation_matrix_optimized(float x0,float x1,float x2, SMatrix33* matrix)
{
	//input: x0,x1, and x2 are random numbers sampled between 0 and 1
	//output: a uniformly sampled random 3D rotation matrix

	//
	// Though I trust the book this is derived from I should do some benchmarking just to verify
	// that my implementation is working correctly.
	// Benchmarking would basically involve generating histograms of x,y,z locations on the sphere and checking they line up
	//
	// This version of the random rotation looks to be working
	// at least the determinant and the transpose of the matrix are the same
	// as expected for a rotation matrix
	//
	
	//cout << " x: " << x0 << " y: " << x1 << " z: " << x2 << endl;

	//Adapted from Graphics Gems III Fast Random Rotation Matrices, Section by James Avro, appendix includes C code algorithm
	Double_t theta,phi,z;
	Double_t r,st,ct;
	Double_t Sx,Sy;
	Double_t Vx, Vy, Vz;

	//Normalize the input random numbers to new values
	theta = x0 * TMath::Pi() * 2;
	phi = x1 * TMath::Pi() * 2;
	z = x2 * 2;

	r = TMath::Sqrt(z);
	Vx = TMath::Sin(phi)*r;
	Vy = TMath::Cos(phi)*r;
	Vz = TMath::Sqrt(2-z);

	st = TMath::Sin(theta);
	ct = TMath::Cos(theta);
	Sx = Vx * ct - Vy * st;
	Sy = Vx * st + Vy * ct;

	double matrixvals[9] = {Vx*Sx-ct,Vx*Sy-st,Vx*Vz,Vy*Sx+st,Vy*Sy-ct,Vy*Vz,Vz*Sx,Vz*Sy,1-z};

	matrix->SetElements(matrixvals,matrixvals+9);
}