#include "Track.cuh"
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TEntryList.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TMath.h"
#include "ThreadTask.hh"
#include <tuple>

Track::Track() {}

void Track::AllocateAndLoadTrack(const ThreadTask& task)
{
	long nVals = task.GetExitPoint() - task.GetEntryPoint() + 1; //+1 because number of values includes first and last value
	size_t trackSize = nVals * sizeof(double);

	cudaMalloc(&x,nVals*sizeof(double));
	cudaMalloc(&y,nVals*sizeof(double));
	cudaMalloc(&z,nVals*sizeof(double));
	cudaMalloc(&edep,nVals*sizeof(double));

	//Open the file in each process and make a Tree Reader
	TFile f = TFile(task.GetFilename());
	TTreeReader trackReader("Tracks", &f);
	trackReader.SetEntriesRange(task.GetEntryPoint(),task.GetExitPoint());
	TTreeReaderValue<double_t> xReader(trackReader, "x [nm]");
	TTreeReaderValue<double_t> yReader(trackReader, "y [nm]");
	TTreeReaderValue<double_t> zReader(trackReader, "z [nm]");
	TTreeReaderValue<double_t> edepReader(trackReader, "edep [eV]");

	std::cout << "thread #: " << task.GetThreadID() << " starting at: " << std::to_string(task.GetEntryPoint()) << std::endl;

	//malloc and cudaMalloc our arrays respectively
	double* xhost = (double *)malloc(trackSize);
	double* yhost = (double *)malloc(trackSize);
	double* zhost = (double *)malloc(trackSize);
	double* edephost = (double *)malloc(trackSize);

	//Fill the unified memory arrays from the CPU
	for (long loopnum = 0; trackReader.Next(); loopnum++) 
	{
		xhost[loopnum] = *xReader;
		yhost[loopnum] = *yReader;
		zhost[loopnum] = *zReader;
		edephost[loopnum] = *edepReader;
	}

	//Copy track to GPU memory
	cudaMemcpy(x,xhost,trackSize,cudaMemcpyHostToDevice);
	cudaMemcpy(y,yhost,trackSize,cudaMemcpyHostToDevice);
	cudaMemcpy(z,zhost,trackSize,cudaMemcpyHostToDevice);
	cudaMemcpy(edep,edephost,trackSize,cudaMemcpyHostToDevice);

	//Free host track
	free(xhost);
	free(yhost);
	free(zhost);
	free(edephost);
}

void Track::AllocateAndLoadTrack(const std::tuple<Int_t,Int_t,Int_t,TString>& input)
{
	long nVals = std::get<1>(input) - std::get<0>(input) + 1; //+1 because number of values includes first and last value
	size_t trackSize = nVals * sizeof(double);

	cudaMalloc(&x,nVals*sizeof(double));
	cudaMalloc(&y,nVals*sizeof(double));
	cudaMalloc(&z,nVals*sizeof(double));
	cudaMalloc(&edep,nVals*sizeof(double));

	//Open the file in each process and make a Tree Reader
	TFile f = TFile(std::get<3>(input));
	TTreeReader trackReader("Tracks", &f);
	trackReader.SetEntriesRange(std::get<0>(input),std::get<1>(input));
	TTreeReaderValue<double_t> xReader(trackReader, "x [nm]");
	TTreeReaderValue<double_t> yReader(trackReader, "y [nm]");
	TTreeReaderValue<double_t> zReader(trackReader, "z [nm]");
	TTreeReaderValue<double_t> edepReader(trackReader, "edep [eV]");

	std::cout << "thread #: " << std::get<2>(input) << " starting at: " << std::to_string(std::get<0>(input)) << std::endl;

	//malloc and cudaMalloc our arrays respectively
	double* xhost = (double *)malloc(trackSize);
	double* yhost = (double *)malloc(trackSize);
	double* zhost = (double *)malloc(trackSize);
	double* edephost = (double *)malloc(trackSize);

	//Fill the unified memory arrays from the CPU
	for (long loopnum = 0; trackReader.Next(); loopnum++) 
	{
		xhost[loopnum] = *xReader;
		yhost[loopnum] = *yReader;
		zhost[loopnum] = *zReader;
		edephost[loopnum] = *edepReader;
	}

	//Copy track to GPU memory
	cudaMemcpy(x,xhost,trackSize,cudaMemcpyHostToDevice);
	cudaMemcpy(y,yhost,trackSize,cudaMemcpyHostToDevice);
	cudaMemcpy(z,zhost,trackSize,cudaMemcpyHostToDevice);
	cudaMemcpy(edep,edephost,trackSize,cudaMemcpyHostToDevice);

	//Free host track
	free(xhost);
	free(yhost);
	free(zhost);
	free(edephost);
}

void Track::AllocateEmptyTrack(int nVals)
{
	cudaMalloc(&x,nVals*sizeof(double));
	cudaMalloc(&y,nVals*sizeof(double));
	cudaMalloc(&z,nVals*sizeof(double));
	cudaMalloc(&edep,nVals*sizeof(double));
}

void Track::Free()
{
	cudaFree(x);
	cudaFree(y);
	cudaFree(z);
	cudaFree(edep);
}

