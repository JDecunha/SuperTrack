//SuperTrack
#include "CPUSuperimpose.hh"
#include "GPUSuperimpose.cuh"
#include "testCUDA.cuh"
#include "utils.hh"
//ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TPad.h"


//SMatrix and SVector are the fastest
//ways to hold vectors and matrices in ROOT
typedef ROOT::Math::SVector<Double_t,3> SVector3;
#define VERBOSE 0
#define VERY_VERBOSE 1

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

	TH1F histo = score_lineal_GPU("/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/software/MicroTrackGenerator/output/proton/50.0MeV/4060394578999227944_thread_0.root",5e3,5e3,2,2000);

	int end_time = time(0); cout << "ROOT Program Ending. Seconds elapsed: " << (end_time-start_time) << endl;	
}

void SuperTrack()
{
	GPU_lineal_test();
}

# ifndef __CINT__
int main()
{
  SuperTrack();
  return 0;
}
# endif


