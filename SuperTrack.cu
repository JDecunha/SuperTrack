//SuperTrack
#include "CPUSuperimpose.hh"
#include "GPUSuperimpose.cuh"
#include "utils.hh"
#include "ThreadAllocator.hh"
#include "VoxelConstrainedSphereMethod.hh"
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
	//Create the macro file
	INIReader reader = INIReader("../macros/test.ini");
	
	//Allocate tasks for threads
	std::string folderPath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/software/MicroTrackGenerator/output/proton/50.0MeV/"; 
	ThreadAllocator folderAllocator = ThreadAllocator(folderPath,4,100,1,2);
	std::vector<ThreadAllocation> ThreadAllocations;	
	folderAllocator.ReturnThreadAllocations(ThreadAllocations);

	//Create the simulation method based off of information in the reader
	SimulationMethod* method = new VoxelConstrainedSphereMethod(reader);

	score_lineal_GPU_New(ThreadAllocations,5e3,5e3);

	/*
	So how will the SuperTrackManager be implemented.

	SuperTrackManager manager = SuperTrackManager();
	manager.AddHistogram(Histogram);
	manager.AddSimulationMethod(method);
	manager.AddThreadAllocations(ThreadAllocations);
	manager.Run();
	*/
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


