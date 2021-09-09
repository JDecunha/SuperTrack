//SuperTrack
#include "testCUDA.cuh"
#include "utils.hh"
//ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TEntryList.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TMath.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TPad.h"
#include "ROOT/TProcessExecutor.hxx"
//STD
#include <unordered_map>
#include <vector>
#include <iterator>
#include <tuple>
#include <filesystem>
//CUDA Libraries
#include <cuda.h>
#include <curand.h>

//SMatrix and SVector are the fastest
//ways to hold vectors and matrices in ROOT
typedef ROOT::Math::SVector<Double_t,3> SVector3;
#define VERBOSE 0
#define VERY_VERBOSE 1

using namespace std;

TH1F score_lineal_GPU(TString filepath, float_t scoring_sphere_spacing, float_t scoring_sphere_diameter, Int_t nthreads, Int_t nSamples = 1, Long_t random_seed = time(NULL))
{
	//open the file and retrieve the trees
	TFile f = TFile(filepath);
	TTree *trackIndex;
	f.GetObject("Track index",trackIndex);
	long long nTracksInFile = trackIndex->GetEntries();

	//Populate our tuple with the first entry, last entry, and random seed for each thread
	std::vector<std::tuple<Int_t,Int_t,Int_t,TString>> perthread_input_arguments;

	//TODO: update the GEANT code to use long long for the event index
	if (nTracksInFile <= nthreads)
	{ 
		long start_entry_val = 0;
		TTreeReader trackIndexReader("Track index", &f);
		TTreeReaderValue<long long> end_entry_val(trackIndexReader, "index");

		for (Int_t i = 0; i < nTracksInFile; i++)
		{
			trackIndexReader.Next();
			perthread_input_arguments.push_back(std::make_tuple(start_entry_val,*end_entry_val-1,i,filepath));
			//Wcout << "thread: " << i << " start val: " << start_entry_val << " end val: " << *end_entry_val-1 << endl;
			start_entry_val = *end_entry_val;
		}
	}
	else
	{
		cout << "Number of tracks in file greater than requested threads. Case not yet implemented." << endl;
	}

	//We are done reading the Tree single threaded. Close it.
	f.Close();


	
	//the = sign captures everything in the enclosing function by value. Meaning it makes a process local copy.
	auto workItem = [=](std::tuple<Int_t,Int_t,Int_t,TString> input) 
	{
		//Open the file in each process and make a Tree Reader
		TFile f = TFile(get<3>(input));
		TTreeReader trackReader("Tracks", &f);
		trackReader.SetEntriesRange(get<0>(input),get<1>(input));
		TTreeReaderValue<double_t> xReader(trackReader, "x [nm]");
		TTreeReaderValue<double_t> yReader(trackReader, "y [nm]");
		TTreeReaderValue<double_t> zReader(trackReader, "z [nm]");
		TTreeReaderValue<double_t> edepReader(trackReader, "edep [eV]");

		cout << "thread #: " << get<2>(input) << " starting at: " << to_string(get<0>(input)) << endl;

		//Determine size of arrays. Define them. Then allcate unified memory on CPU and GPU
		long nVals = get<1>(input) - get<0>(input) + 1; //+1 because number of values includes first and last value
		size_t trackSize = nVals * sizeof(double);

		double *x;
		double *y;
		double *z;
		double *edep;

		cudaMallocManaged(&x,trackSize);
		cudaMallocManaged(&y,trackSize);
		cudaMallocManaged(&z,trackSize);
		cudaMallocManaged(&edep,trackSize);

		//Fill the unified memory arrays from the CPU
		for (long loopnum = 0; trackReader.Next(); loopnum++) 
		{
			x[loopnum] = *xReader;
			y[loopnum] = *yReader;
			z[loopnum] = *zReader;
			edep[loopnum] = *edepReader;
		}

		//Get the device Id for active GPU
		int deviceId;
		cudaGetDevice(&deviceId);    
		//Prefetch memory by the GPU
		cudaMemPrefetchAsync(x,trackSize,deviceId);
		cudaMemPrefetchAsync(y,trackSize,deviceId);
		cudaMemPrefetchAsync(z,trackSize,deviceId);
		cudaMemPrefetchAsync(edep,trackSize,deviceId);

		//Allocate GPU only memory for the random numbers
		float *randomVals;
		cudaMalloc(&randomVals,2*sizeof(float)*nSamples); //2 values for x,y times the number of oversamples needed
		
		//Random number generation on GPU
		curandGenerator_t randGenerator;
		curandCreateGenerator(&randGenerator,CURAND_RNG_PSEUDO_DEFAULT); //consider changing this to Mersenne Twister later
		curandSetPseudoRandomGeneratorSeed(randGenerator,random_seed+get<2>(input));
		curandGenerateUniform(randGenerator,randomVals,2*nSamples);

		//Allocate GPU only memory for the volume:edep paired list
		long *volumeID;
		double *edepInVolume;

		cudaMalloc(&volumeID,sizeof(long)*nVals);
		cudaMalloc(&edepInVolume,trackSize);

		//TODO: Invoke superimposing kernel call here

		//TODO: Consolidate results of superimposing into a hash table

		//TODO: Transform hash table into a histogram

		//TODO: Transfer histogram back to CPU memory and return


	  	//Initialize the histogram
		TH1F lineal_histogram = TH1F("Lineal energy histogram", "y*f(y)", 200, -2,1);
		return lineal_histogram;

	};


	// Create the pool of workers
  ROOT::TProcessExecutor workers(nthreads);
  //Process the jobs and get a vector of the output
  std::vector<TH1F> process_output = workers.Map(workItem, perthread_input_arguments);

   TH1F lineal_histogram = TH1F("Lineal energy histogram", "y*f(y)", 200, -2,1);

   return lineal_histogram;

}

