//SuperTrack
#include "CPUSuperimpose.hh"
#include "GPUSuperimpose.cuh"
#include "utils.hh"
#include "ThreadAllocator.hh"
#include "VoxelConstrainedSphereMethod.hh"
#include "SuperTrackManager.hh"
#include "SimulationMethodFactory.hh"
//inih
#include "INIReader.h"
//ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TPad.h"

using namespace std;


void CPU_lineal_test()
{
	int start_time = time(0); cout << "ROOT Program Beginning" << endl; 	

	TH1F histo = score_lineal_voxel("/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/software/MicroTrackGenerator/output/proton/50.0MeV/4060394578999227944_thread_0.root",5e3,5e3,2,50);

	//Plotting
	TCanvas *c = new TCanvas();
	THStack *histo_stack = new THStack("histograms","");
	histo_stack->Add(&histo);
	histo_stack->Draw("nostack"); //Draw histogram
	gPad->SetLogx(); //Set the logarithmic axes appropriately
	gPad->Modified(); 
	c->Print("mostrecenthist.png");

	int end_time = time(0); cout << "ROOT Program Ending. Seconds elapsed: " << (end_time-start_time) << endl;	
}

void GPU_lineal_test()
{
	int start_time = time(0); cout << "ROOT Program Beginning" << endl; 	

	TH1F histo = score_lineal_GPU("/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/software/MicroTrackGenerator/output/proton/50.0MeV/4060394578999227944_thread_0.root",5e3,5e3,2,100);

	int end_time = time(0); cout << "ROOT Program Ending. Seconds elapsed: " << (end_time-start_time) << endl;	
}

void File_Allocator_test()
{
	//Read the macro file
	INIReader reader = INIReader("../macros/test.ini");

	//Add the voxel constrained sphere method to the list of available simulation methods
	SimulationMethodFactory& methodFactory = SimulationMethodFactory::GetInstance();

	methodFactory.AddSimulationMethod("VoxelConstrainedSphere", &VoxelConstrainedSphereMethod::Construct);

	//Construct the simulation method based upon the information in the macro file
	SimulationMethod* method = methodFactory.Construct(reader);

	//Let's think about the design philosophy though
	//1.) SimulationMethodFactory: is a thread-safe singleton. Eagerly initialized. Different threads will access it but no race conditions should arise since it's just giving you a pointer to a constructor.
	//2.) SimulationMethod: Should be initialized once on each thread, since it holds memory allocations that need to be distinct on each thread.
	//3.) Histogram: Should be initialized once on each thread, for the same reason as SimulationMethod.
	//4.) ThreadAllocations: Should be initiated on the main thread, and sent to SuperTrackManager.
	//5.) SuperTrackManager: Should be another singleton. It takes the .ini path as an input. Creates the Thread Allocations. Then creates a simulation method and histogram on every single thread.

	/*
	//Allocate tasks for threads
	std::string folderPath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/software/MicroTrackGenerator/output/proton/50.0MeV/"; 
	ThreadAllocator folderAllocator = ThreadAllocator(folderPath,4,100,1,2);
	std::vector<ThreadAllocation> ThreadAllocations;	
	folderAllocator.ReturnThreadAllocations(ThreadAllocations);

	//Make histogram
	Histogram* histogram = new Histogram(200,-1,2,"log");

	//Create the simulation method based off of information in the reader
	SimulationMethod* method = new VoxelConstrainedSphereMethod(reader);

	//score_lineal_GPU_New(ThreadAllocations,5e3,5e3);

	SuperTrackManager manager = SuperTrackManager();
	manager.AddThreadAllocations(ThreadAllocations);
	manager.AddHistogram(histogram);
	manager.AddSimulationMethod(method);

	manager.Initialize();
	manager.Run();*/
	
}

void SuperTrack()
{
	/*INIReader reader("../macros/test.ini");

	if (reader.ParseError() < 0) 
	    std::cout << "Can't load 'test.ini'\n";
	
	cout << reader.Get("user", "name", "UNKNOWN") << endl;*/

	File_Allocator_test();
}

# ifndef __CINT__
int main()
{
  SuperTrack();
  return 0;
}
# endif


