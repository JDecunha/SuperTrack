#include "Track.cuh"
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TEntryList.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TMath.h"
#include <tuple>

//TODO: this is a bad design principle
//replace this with separate allocation function
Track::Track(int nVals)
{
	cudaMalloc(&x,nVals*sizeof(double));
	cudaMalloc(&y,nVals*sizeof(double));
	cudaMalloc(&z,nVals*sizeof(double));
	cudaMalloc(&edep,nVals*sizeof(double));
}

Track::~Track()
{
	cudaFree(&x);
	cudaFree(&y);
	cudaFree(&z);
	cudaFree(&edep);
}

void Track::LoadTrack(const std::tuple<Int_t,Int_t,Int_t,TString> &input, Track *deviceTrack)
{
	//Open the file in each process and make a Tree Reader
	TFile f = TFile(std::get<3>(input));
	TTreeReader trackReader("Tracks", &f);
	trackReader.SetEntriesRange(std::get<0>(input),std::get<1>(input));
	TTreeReaderValue<double_t> xReader(trackReader, "x [nm]");
	TTreeReaderValue<double_t> yReader(trackReader, "y [nm]");
	TTreeReaderValue<double_t> zReader(trackReader, "z [nm]");
	TTreeReaderValue<double_t> edepReader(trackReader, "edep [eV]");

	std::cout << "thread #: " << std::get<2>(input) << " starting at: " << std::to_string(std::get<0>(input)) << std::endl;

	//Determine size of arrays. Define them. Then allcate unified memory on CPU and GPU
	long nVals = std::get<1>(input) - std::get<0>(input) + 1; //+1 because number of values includes first and last value
	size_t trackSize = nVals * sizeof(double);

	//malloc and cudaMalloc our arrays respectively
	double* x = (double *)malloc(trackSize);
	double* y = (double *)malloc(trackSize);
	double* z = (double *)malloc(trackSize);
	double* edep = (double *)malloc(trackSize);

	//Fill the unified memory arrays from the CPU
	for (long loopnum = 0; trackReader.Next(); loopnum++) 
	{
		x[loopnum] = *xReader;
		y[loopnum] = *yReader;
		z[loopnum] = *zReader;
		edep[loopnum] = *edepReader;
	}

	//Copy track to GPU memory
	cudaMemcpy(deviceTrack->x,x,trackSize,cudaMemcpyHostToDevice);
	cudaMemcpy(deviceTrack->y,y,trackSize,cudaMemcpyHostToDevice);
	cudaMemcpy(deviceTrack->z,z,trackSize,cudaMemcpyHostToDevice);
	cudaMemcpy(deviceTrack->edep,edep,trackSize,cudaMemcpyHostToDevice);

	//TODO: free host track, maybe we can't free it yet because it's still being copied right
}

