//SuperTrack
#include "SuperTrackManager.hh"
#include "VoxelConstrainedSphereMethod.hh"
//inih
#include "INIReader.h"
//gTest
#include "gtest/gtest.h"
//std
#include <iostream>

TEST(SuperTrackVoxelConstrainedTest, A1) //test 1 edep, in 1 volume, at 1 keV/um
{
    //Step 1.) Set the macro filepath for this test
    std::string INIPath = "./testMacros/testA1.ini";

    //Step 2.) Run SuperTrack
    SuperTrackManager& manager = SuperTrackManager::GetInstance();
    manager.AddSimulationMethod("VoxelConstrainedSphere", &VoxelConstrainedSphereMethod::Construct);
    manager.Initialize(new INIReader(INIPath));
    manager.Run();
    TH1D output = manager.GetOutput();

    ///
    /// Step 3.) Test the output
    ///

    double edepValue = 1; //keV/um
    int numEdeps = 1; //1 energy deposition event

    int nBins = output.GetNbinsX(); //get histogram nbins

    for (int i = 1; i <= nBins; i++) //Loop over each bin
    {
        auto counts = output.GetBinContent(i);
        auto low_edge = output.GetBinLowEdge(i);
        auto high_edge = output.GetBinLowEdge(i+1);
        
        //Counts should equal numEdeps in the specified bin
        if(low_edge <= edepValue && high_edge > edepValue) //<= because lower bin is inclusive
        { 
            EXPECT_EQ(counts, numEdeps) << "Counts equal to " << counts << " rather than " << numEdeps << " in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        } 
        else //Counts should be 0 in all bins except the specified bin
        {
            EXPECT_EQ(counts,0) << "Counts equal to " << counts << " rather than 0 in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        }
    }
}

TEST(SuperTrackVoxelConstrainedTest, A2) //test 1 edep, in every volume in a line, at 1 keV/um
{
    //Step 1.) Set the macro filepath for this test
    std::string INIPath = "./testMacros/testA2.ini";

    //Step 2.) Run SuperTrack
    SuperTrackManager& manager = SuperTrackManager::GetInstance();
    manager.AddSimulationMethod("VoxelConstrainedSphere", &VoxelConstrainedSphereMethod::Construct);
    manager.Initialize(new INIReader(INIPath));
    manager.Run();
    TH1D output = manager.GetOutput();

    ///
    /// Step 3.) Test the output
    ///

    double edepValue = 1; //keV/um
    int numEdeps = 3000; //Every single sphere in a plane should have an edep in it (3000^2)

    int nBins = output.GetNbinsX(); //get histogram nbins

    for (int i = 1; i <= nBins; i++) //Loop over each bin
    {
        auto counts = output.GetBinContent(i);
        auto low_edge = output.GetBinLowEdge(i);
        auto high_edge = output.GetBinLowEdge(i+1);
        
        //Counts should equal numEdeps in the specified bin
        if(low_edge <= edepValue && high_edge > edepValue) //<= because lower bin is inclusive
        { 
            EXPECT_EQ(counts, numEdeps) << "Counts equal to " << counts << " rather than " << numEdeps << " in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        } 
        else //Counts should be 0 in all bins except the specified bin
        {
            EXPECT_EQ(counts,0) << "Counts equal to " << counts << " rather than 0 in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        }
    }
}

TEST(SuperTrackVoxelConstrainedTest, A3) //test 1 edep, in every volume in a plane of the cube, at 1 keV/um
{
    //Step 1.) Set the macro filepath for this test
    std::string INIPath = "./testMacros/testA3.ini";

    //Step 2.) Run SuperTrack
    SuperTrackManager& manager = SuperTrackManager::GetInstance();
    manager.AddSimulationMethod("VoxelConstrainedSphere", &VoxelConstrainedSphereMethod::Construct);
    manager.Initialize(new INIReader(INIPath));
    manager.Run();
    TH1D output = manager.GetOutput();

    ///
    /// Step 3.) Test the output
    ///

    double edepValue = 1; //keV/um
    int numEdeps = 9000000; //Every single sphere in a plane should have an edep in it (3000^2)

    int nBins = output.GetNbinsX(); //get histogram nbins

    for (int i = 1; i <= nBins; i++) //Loop over each bin
    {
        auto counts = output.GetBinContent(i);
        auto low_edge = output.GetBinLowEdge(i);
        auto high_edge = output.GetBinLowEdge(i+1);
        
        //Counts should equal numEdeps in the specified bin
        if(low_edge <= edepValue && high_edge > edepValue) //<= because lower bin is inclusive
        { 
            EXPECT_EQ(counts, numEdeps) << "Counts equal to " << counts << " rather than " << numEdeps << " in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        } 
        else //Counts should be 0 in all bins except the specified bin
        {
            EXPECT_EQ(counts,0) << "Counts equal to " << counts << " rather than 0 in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        }
    }
}

TEST(SuperTrackVoxelConstrainedTest, A4) //same as test A2 but with 5 tracks and 1000 oversamples per track
{
    //Step 1.) Set the macro filepath for this test
    std::string INIPath = "./testMacros/testA4.ini";

    //Step 2.) Run SuperTrack
    SuperTrackManager& manager = SuperTrackManager::GetInstance();
    manager.AddSimulationMethod("VoxelConstrainedSphere", &VoxelConstrainedSphereMethod::Construct);
    manager.Initialize(new INIReader(INIPath));
    manager.Run();
    TH1D output = manager.GetOutput();

    ///
    /// Step 3.) Test the output
    ///

    double edepValue = 1; //keV/um
    int numEdeps = 3000*1000*5; //Every single sphere in a plane should have an edep in it (3000^2)

    int nBins = output.GetNbinsX(); //get histogram nbins

    for (int i = 1; i <= nBins; i++) //Loop over each bin
    {
        auto counts = output.GetBinContent(i);
        auto low_edge = output.GetBinLowEdge(i);
        auto high_edge = output.GetBinLowEdge(i+1);
        
        //Counts should equal numEdeps in the specified bin
        if(low_edge <= edepValue && high_edge > edepValue) //<= because lower bin is inclusive
        { 
            EXPECT_EQ(counts, numEdeps) << "Counts equal to " << counts << " rather than " << numEdeps << " in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        } 
        else //Counts should be 0 in all bins except the specified bin
        {
            EXPECT_EQ(counts,0) << "Counts equal to " << counts << " rather than 0 in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        }
    }
}

TEST(SuperTrackVoxelConstrainedTest, A5) //5 tracks, with 1000 edeps in each volume, oversampled 1000 times
{
    //A5: Have a test that tests multiple energy depositions from one track within a single sphere
    //A5: Have a test that tests the sensitivity of spheres to 'see' energy depositions within them (From the edge)
    //For this test we have a safety parameter of zero! and we still see the edeps

    //Step 1.) Set the macro filepath for this test
    std::string INIPath = "./testMacros/testA5.ini";

    //Step 2.) Run SuperTrack
    SuperTrackManager& manager = SuperTrackManager::GetInstance();
    manager.AddSimulationMethod("VoxelConstrainedSphere", &VoxelConstrainedSphereMethod::Construct);
    manager.Initialize(new INIReader(INIPath));
    manager.Run();
    TH1D output = manager.GetOutput();

    ///
    /// Step 3.) Test the output
    ///

    double edepValue = 1000; //keV/um
    int numEdeps = 5*1000*3000; //3000 spheres in a row, 1000 oversamples, for 5 distinct tracks

    int nBins = output.GetNbinsX(); //get histogram nbins

    for (int i = 1; i <= nBins; i++) //Loop over each bin
    {
        auto counts = output.GetBinContent(i);
        auto low_edge = output.GetBinLowEdge(i);
        auto high_edge = output.GetBinLowEdge(i+1);
        
        //Counts should equal numEdeps in the specified bin
        if(low_edge <= edepValue && high_edge > edepValue) //<= because lower bin is inclusive
        { 
            EXPECT_EQ(counts, numEdeps) << "Counts equal to " << counts << " rather than " << numEdeps << " in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        } 
        else //Counts should be 0 in all bins except the specified bin
        {
            EXPECT_EQ(counts,0) << "Counts equal to " << counts << " rather than 0 in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        }
    }
}

TEST(SuperTrackVoxelConstrainedTest, A6) //All the edeps are on the surface of the cube
{
    //A6: A test that places edeps on the surface of the box and checks that none of them are detected

    //Step 1.) Set the macro filepath for this test
    std::string INIPath = "./testMacros/testA6.ini";

    //Step 2.) Run SuperTrack
    SuperTrackManager& manager = SuperTrackManager::GetInstance();
    manager.AddSimulationMethod("VoxelConstrainedSphere", &VoxelConstrainedSphereMethod::Construct);
    manager.Initialize(new INIReader(INIPath));
    manager.Run();
    TH1D output = manager.GetOutput();

    ///
    /// Step 3.) Test the output
    ///

    int nBins = output.GetNbinsX(); //get histogram nbins

    for (int i = 1; i <= nBins; i++) //Loop over each bin. There should be zero edeps anywhere
    {
        auto counts = output.GetBinContent(i);
        auto low_edge = output.GetBinLowEdge(i);
        auto high_edge = output.GetBinLowEdge(i+1);
        
        EXPECT_EQ(counts,0) << "Counts equal to " << counts << " rather than 0 in bin of low edge: " << low_edge << " and high edge: " << high_edge;
    }
}

TEST(SuperTrackVoxelConstrainedTest, A7) 
{

    //A test such that you can calculate pi from the randomness of your tracks

    //Step 1.) Set the macro filepath for this test
    std::string INIPath = "./testMacros/testA7.ini";

    //Step 2.) Run SuperTrack
    SuperTrackManager& manager = SuperTrackManager::GetInstance();
    manager.AddSimulationMethod("VoxelConstrainedSphere", &VoxelConstrainedSphereMethod::Construct);
    manager.Initialize(new INIReader(INIPath));
    manager.Run();
    TH1D output = manager.GetOutput();

    ///
    /// Step 3.) Test the output
    ///

    long nTotal = long(20000)*long(10000); //20000 edeps total, 1000 oversamples, 1 track 
    double nEdepsInSpheres = 0;

    int nBins = output.GetNbinsX(); //get histogram nbins

    for (int i = 1; i <= nBins; i++) //Loop over each bin
    {
        long counts = output.GetBinContent(i);
        double low_edge = output.GetBinLowEdge(i);
        double high_edge = output.GetBinLowEdge(i+1);
        double mid_bin = (high_edge+low_edge)/2.;

        double nEdeps = mid_bin*counts;
        
        nEdepsInSpheres += nEdeps; 
    }

    double pi_estimate = nEdepsInSpheres*6/nTotal;
    EXPECT_TRUE(pi_estimate < 3.3 && pi_estimate > 3.05) << "Pi estimate equal to: " << pi_estimate << ". Try running again. This test has a finite probability to fail depending on the random seed.";
}

TEST(SuperTrackVoxelConstrainedTest, A8) //test whether overflow occurs at 4 byte int value
{
    //Step 1.) Set the macro filepath for this test
    std::string INIPath = "./testMacros/testA8.ini";

    //Step 2.) Run SuperTrack
    SuperTrackManager& manager = SuperTrackManager::GetInstance();
    manager.AddSimulationMethod("VoxelConstrainedSphere", &VoxelConstrainedSphereMethod::Construct);
    manager.Initialize(new INIReader(INIPath));
    manager.Run();
    TH1D output = manager.GetOutput();

    ///
    /// Step 3.) Test the output
    ///

    double edepValue = 1; //keV/um
    int numEdeps = 7; //1 energy deposition event

    int nBins = output.GetNbinsX(); //get histogram nbins

    for (int i = 1; i <= nBins; i++) //Loop over each bin
    {
        auto counts = output.GetBinContent(i);
        auto low_edge = output.GetBinLowEdge(i);
        auto high_edge = output.GetBinLowEdge(i+1);
        
        //Counts should equal numEdeps in the specified bin
        if(low_edge <= edepValue && high_edge > edepValue) //<= because lower bin is inclusive
        { 
            EXPECT_EQ(counts, numEdeps) << "Counts equal to " << counts << " rather than " << numEdeps << " in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        } 
        else //Counts should be 0 in all bins except the specified bin
        {
            EXPECT_EQ(counts,0) << "Counts equal to " << counts << " rather than 0 in bin of low edge: " << low_edge << " and high edge: " << high_edge;
        }
    }
}

