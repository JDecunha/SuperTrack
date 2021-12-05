//SuperTrack
#include "SuperTrackManager.hh"
#include "VoxelConstrainedSphereMethod.hh"
//inih
#include "INIReader.h"

void SuperTrack()
{
	time_t start;
	time_t end;

	time(&start);

	//Get the SuperTrackManager
	SuperTrackManager& manager = SuperTrackManager::GetInstance();

	//Add the voxel constrained sphere method to the available simulation methods
	manager.AddSimulationMethod("VoxelConstrainedSphere", &VoxelConstrainedSphereMethod::Construct);

	//Initialize the SuperTrackManager with the input file
	manager.Initialize(new INIReader("../macros/test.ini"));
	manager.Run();

	time(&end);
	std::cout << "SuperTrack concluding after: " << end-start << " seconds." << std::endl;
}

# ifndef __CINT__
int main()
{
  SuperTrack();
  return 0;
}
# endif


//Design philosophy
//1.) SuperTrackManager: Is a singleton. It takes the .ini path as an input. Creates the Thread Allocations. Then creates a simulation method and histogram on every single thread.
//2.) ThreadAllocations: Will be initiated on the main thread and belong to SuperTrackManager.
//3.) SimulationMethodFactory: Is a singleton. Eagerly initialized. Different threads will access it but no race conditions should arise since it's just giving you a pointer to a constructor.
//4.) SimulationMethod: Will be initialized once on each thread, since it holds memory allocations that need to be distinct on each thread.
//5.) Histogram: Will be initialized once on each thread, since it holds memory allocations that need to be distinct on each thread.

//Plotting
/*
TCanvas *c = new TCanvas();
THStack *histo_stack = new THStack("histograms","");
histo_stack->Add(&output_reduced);
histo_stack->Draw("nostack"); //Draw histogram
gPad->SetLogx(); //Set the logarithmic axes appropriately
gPad->Modified(); 
c->Print("mostrecenthist.png");*/