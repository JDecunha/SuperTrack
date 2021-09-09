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
		std::cout << "Number of tracks in file greater than requested threads. Case not yet implemented." << std::endl;
	}

	//We are done reading the Tree single threaded. Close it.
	f.Close();


	
	//the = sign captures everything in the enclosing function by value. Meaning it makes a process local copy.
	auto workItem = [=](std::tuple<Int_t,Int_t,Int_t,TString> input) 
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
		curandSetPseudoRandomGeneratorSeed(randGenerator,random_seed+std::get<2>(input));
		curandGenerateUniform(randGenerator,randomVals,2*nSamples);
		curandDestroyGenerator(randGenerator);

		//Allocate GPU only memory for the volume:edep paired list
		long *volumeID;
		double *edepInVolume;

		cudaMalloc(&volumeID,sizeof(long)*nVals);
		cudaMalloc(&edepInVolume,trackSize);

		//TODO: Invoke superimposing kernel call here

		//TODO: Consolidate results of superimposing into a hash table

		//TODO: Transform hash table into a histogram

		//TODO: Transfer histogram back to CPU memory and return

		//


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

//TODO: Change this to work with a C-style struct later, so x,y,z,edep are all one entry
__global__ void SuperimposeTrack(double greatestSphereOffset, double sphereDiameter, long numSpheresLinear, double* randomVals, double* x, double* y, double* z, double* edep,long *volumeID, double *edepOutput, long numElements)
{
	//Our entire geometry should be able to be described by only the greatest offset, the sphere diameter and number of spheres in a line. That's useful
	double sphereRadius = sphereDiameter/2;
	double linealDenominator = (2./3.)*sphereDiameter; //calculate this here as an efficiency gain

	//Determine index and stride
 	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;

	//Loop over all the energy deposition points
	for (long i = index; i < numElements; i+=stride)
	{
		//Write a zero to edepOutput and volumeID. Doing this here avoids warp divergence later.
		edepOutput[i] = 0; volumeID[i] = 0;

		//Apply random shift
		double x_shifted = x[i] + randomVals[(i*2)];
		double y_shifted = y[i] + randomVals[(i*2)+1];

		//Determine inside box
		if (abs(x_shifted) < abs(greatestSphereOffset)+(sphereRadius) && abs(y_shifted) < abs(greatestSphereOffset)+(sphereRadius) && abs(z[i]) < abs(greatestSphereOffset)+(sphereRadius))
		{
			//Convert position to index in the grid of spheres
			long xIndex = llround((x_shifted-greatestSphereOffset)/sphereDiameter);
			long yIndex = llround((y_shifted-greatestSphereOffset)/sphereDiameter);
			long zIndex = llround((z[i]-greatestSphereOffset)/sphereDiameter);
			
			//Determine the location of the nearest sphere in the grid (with 0,0,0 being the top left sphere, different coordinate system than the ptcls are in)
			double nearestSphereX = xIndex*sphereDiameter;
			double nearestSphereY = yIndex*sphereDiameter;
			double nearestSphereZ = zIndex*sphereDiameter;

			//Find the distance from the nearest sphere. You have to shift x_shift by gSO to get in the same coordinate system as the sphere grid
			//An aside: I feel like there is probably a way that you could define the sphere grid that might reduce the complexity of this kernel
			//Another aside: calculating in cubes would reduce complexity as well
			double distFromNearestSphereX = nearestSphereX-(x_shifted-greatestSphereOffset);
			double distFromNearestSphereY = nearestSphereY-(y_shifted-greatestSphereOffset); 
			double distFromNearestSphereZ = nearestSphereZ-(z[i]-greatestSphereOffset); 

			//Determine if inside the nearest sphere
			double dist = pow(distFromNearestSphereX,2)+pow(distFromNearestSphereY,2)+pow(distFromNearestSphereZ,2);
			dist = sqrt(dist);

			if (dist <= sphereRadius)
			{
				//Determine the Index of the sphere hit
				long sphereHitIndex = xIndex + yIndex*(numSpheresLinear) + zIndex*pow(numSpheresLinear,2); //Keep in mind that for the index it starts counting at zero

				//Write to volumeID and edepOutput
				volumeID[i] = sphereHitIndex;
				edepOutput[i] = edep[i]/linealDenominator; //this should be ev/nm which is same a kev/um
			}
		}
	}


}