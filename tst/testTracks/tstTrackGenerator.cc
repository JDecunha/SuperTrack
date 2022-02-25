#include <iostream>
#include <math.h>

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"

void tstTrackGenerator()
{
	//Need to establish my variables
	double x,y,z,edep;
	long indexEntry = 0;

	//Create TFile
	auto pTrackOutputFile = new TFile("tstTrack1ev1um.root","RECREATE");

	//Code dragged directly out of MicroTrackGenerator
    auto pTrackOutputTree = new TTree("Tracks","Track information data");

    //Configure the branches
    pTrackOutputTree->Branch("x [nm]",&x,"x/D");
    pTrackOutputTree->Branch("y [nm]",&y,"y/D");
    pTrackOutputTree->Branch("z [nm]",&z,"z/D");
    pTrackOutputTree->Branch("edep [eV]",&edep,"edep/D");

    //Create the Tree
    auto pEventIndexTree = new TTree("Track index","Entry number for the end of each track");

    //Configure the branch
    pEventIndexTree->Branch("index",&indexEntry,"index/L");
    
    //Set edep to known value
    edep = 5000; //To give 1 keV per micrometer in a 1 um diameter sphere

    //X and Y along the central axis of the box
    x = 0; y = 0;

   for (int i = 0; i < 1; i++)
    {
        //Place a series of edeps along the z axis, each jumping 1 um 
        //Start at the edge of the box plus half a diameter length
        //Jump one diameter length every iteration
        for (double outputValue = -1.5e6; outputValue < -1.5e6+0.1e3; outputValue+=0.1e3)
        {
        	z = outputValue;
            pTrackOutputTree->Fill();
            pTrackOutputTree->Fill();
            pTrackOutputTree->Fill();
            pTrackOutputTree->Fill();
        	indexEntry += 4;
        }

        pEventIndexTree->Fill();
    }



    pTrackOutputFile->Write(0,TObject::kWriteDelete);
    pTrackOutputFile->Close();
}