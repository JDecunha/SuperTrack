//SuperTrack
#include "SuperTrackManager.hh"
#include "VoxelConstrainedSphereMethod.hh"
//inih
#include "INIReader.h"
//gTest
#include "gtest/gtest.h"

TEST(SuperTrackVoxelConstrainedTest, A1)
{

    //Step 1.) Set the macro filepath for this test
    std::string INIPath = "../tst/testMacros/testA1.ini";

    //Step 2.) Run SuperTrack
    SuperTrackManager& manager = SuperTrackManager::GetInstance();
    manager.AddSimulationMethod("VoxelConstrainedSphere", &VoxelConstrainedSphereMethod::Construct);
    manager.Initialize(new INIReader(INIPath));
    manager.Run();

    //TString wamp;
    //wamp = "yolo";
    //std::cout<<"just yoloing"<<std::endl;

    //Step 3.) Retrieve the output file and verify validity
    /*TFile* f = new TFile("/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/software/SuperTrack/output/test/test_1ev1645576079.root");

    TH1D* h = (TH1D*)f->Get("Lineal energy histogram");

    int nBins = output.GetNbinsX(); //get histogram nbins

    for (int i = 1; i <= nBins; i++) //Loop over each bin
    {
        auto counts = output.GetBinContent(i);
        auto low_edge = output.GetBinLowEdge(i);
        auto high_edge = output.GetBinLowEdge(i+1);
        
        //Counts should equal numEdeps in the specified bin
        if(low_edge <= edepValue && high_edge > edepValue) //<= because lower bin is inclusive
        { 
            EXPECT_EQ(counts, numEdeps);
        } 
        else //Counts should be 0 in all bins except the specified bin
        {
            EXPECT_EQ(counts,0) << "Counts equal to " << counts << " rather than 0 in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        }

    }


    f->Close();*/
}
