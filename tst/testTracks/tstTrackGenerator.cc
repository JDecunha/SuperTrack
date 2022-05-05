//ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
//std library
#include <iostream>
#include <math.h>

namespace generateTestTracks
{ 
    void GenerateTestA1();
    void GenerateTestA2();
    void GenerateTestA3();
    void GenerateTestA4();
    void GenerateTestA5();
    void GenerateTestA6();
    void GenerateTestA7();
    void GenerateTestA8();
};

void tstTrackGenerator()
{
    generateTestTracks::GenerateTestA1();
    generateTestTracks::GenerateTestA2();
    generateTestTracks::GenerateTestA3();
    generateTestTracks::GenerateTestA4();
    generateTestTracks::GenerateTestA5();
    generateTestTracks::GenerateTestA6();
    generateTestTracks::GenerateTestA7();
    generateTestTracks::GenerateTestA8();
}

void generateTestTracks::GenerateTestA1()
{
    //Need to establish my variables
    double x,y,z,edep;
    long indexEntry = 0;

    //Set track information to known fixed values
    x = 0.5e3; y = 0.5e3; //X and Y along the axis of spheres nearest the center of the box
    edep = 2000/3; //2000 eV to give 1 keV per micrometer in a 1 um diameter sphere

    //Create TFile
    auto pTrackOutputFile = new TFile("./testA1/testA1.root","RECREATE");

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
    
    double startPosition = -1.5e6+0.5e3;
    double stopPosition = startPosition+1e3;
    double positionStride = 1e3;

   for (int i = 0; i < 1; i++)
    {
        //Place a series of edeps along the z axis, each jumping 1 um 
        //Start at the edge of the box plus half a diameter length
        //Jump one diameter length every iteration
        for (double outputValue = startPosition; outputValue < stopPosition; outputValue+=positionStride)
        {
            z = outputValue;
            pTrackOutputTree->Fill();
            indexEntry += 1;
        }

        pEventIndexTree->Fill();
    }

    std::cout << indexEntry << std::endl;

    pTrackOutputFile->Write(0,TObject::kWriteDelete);
    pTrackOutputFile->Close();
}

void generateTestTracks::GenerateTestA2()
{
    //Need to establish my variables
    double x,y,z,edep;
    long indexEntry = 0;

    //Set track information to known fixed values
    x = 0.5e3; y = 0.5e3; //X and Y along the axis of spheres nearest the center of the box
    edep = double(2000)/3; //2000 eV to give 1 keV per micrometer in a 1 um diameter sphere

    //Create TFile
    auto pTrackOutputFile = new TFile("./testA2/testA2.root","RECREATE");

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
    
    double startPosition = -1.5e6+0.5e3;
    double stopPosition = 1.5e6-0.5e3;
    double positionStride = 1e3;

   for (int i = 0; i < 1; i++)
    {
        //Place a series of edeps along the z axis, each jumping 1 um 
        //Start at the edge of the box plus half a diameter length
        //Jump one diameter length every iteration
        for (double outputValue = startPosition; outputValue <= stopPosition; outputValue+=positionStride)
        {
            z = outputValue;
            pTrackOutputTree->Fill();
            indexEntry += 1;
        }

        pEventIndexTree->Fill();
    }

    std::cout << indexEntry << std::endl;

    pTrackOutputFile->Write(0,TObject::kWriteDelete);
    pTrackOutputFile->Close();
}

void generateTestTracks::GenerateTestA3() // Fill a plane of the cube (y, and z vary, X is constant)
{
    //Need to establish my variables
    double x,y,z,edep;
    long indexEntry = 0;

    //Set track information to known fixed values
    x = 0.5e3; //X along the axis of spheres nearest the center of the box
    edep = double(2000)/3; //2000 eV to give 1 keV per micrometer in a 1 um diameter sphere

    //Create TFile
    auto pTrackOutputFile = new TFile("./testA3/testA3.root","RECREATE");

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
    
    double startPosition = -1.5e6+0.5e3;
    double stopPosition = 1.5e6-0.5e3;
    double positionStride = 1e3;

    //Place a plane in y and z
    for (int i = 0; i < 1; i++)
    {
        for (double yCoordinate = startPosition; yCoordinate <= stopPosition; yCoordinate+=positionStride)
        {
            for (double zCoordinate = startPosition; zCoordinate <= stopPosition; zCoordinate+=positionStride)
            {
                y = yCoordinate;
                z = zCoordinate;
                pTrackOutputTree->Fill();
                indexEntry += 1;
            }
        }
        pEventIndexTree->Fill();
    }

    std::cout << indexEntry << std::endl;

    pTrackOutputFile->Write(0,TObject::kWriteDelete);
    pTrackOutputFile->Close();
}

void generateTestTracks::GenerateTestA4() // Generate a track to test that filtering outside the box works
{
    //Need to establish my variables
    double x,y,z,edep;
    long indexEntry = 0;

    //Set track information to known fixed values
    x = 0.5e3; y = 0.5e3; //X and Y along the axis of spheres nearest the center of the box
    edep = double(2000)/3.; //2000 eV to give 1 keV per micrometer in a 1 um diameter sphere

    //Create TFile
    auto pTrackOutputFile = new TFile("./testA4/testA4.root","RECREATE");

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
    
    double startPosition = -1.5e6+0.5e3;
    double stopPosition = 1.5e6-0.5e3;
    double positionStride = 1e3;

   for (int i = 0; i < 5; i++)
    {
        //Place a series of edeps along the z axis, each jumping 1 um 
        //Start at the edge of the box plus half a diameter length
        //Jump one diameter length every iteration
        for (double outputValue = startPosition; outputValue <= stopPosition; outputValue+=positionStride)
        {
            z = outputValue;
            pTrackOutputTree->Fill();
            indexEntry += 1;
        }

        pEventIndexTree->Fill();
    }

    std::cout << indexEntry << std::endl;

    pTrackOutputFile->Write(0,TObject::kWriteDelete);
    pTrackOutputFile->Close();
}

void generateTestTracks::GenerateTestA5() // Generate a track to test that filtering outside the box works
{
    //Need to establish my variables
    double x,y,z,edep;
    long indexEntry = 0;

    //Set track information to known fixed values
    x = 0.5e3; y = 0.5e3; //X and Y along the axis of spheres nearest the center of the box
    edep = double(2000)/double(3); //2000 eV to give 1 keV per micrometer in a 1 um diameter sphere

    //Create TFile
    auto pTrackOutputFile = new TFile("./testA5/testA5.root","RECREATE");

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
    
    double startPosition = -1.5e6+0.5e3;
    double stopPosition = 1.5e6-0.5e3;
    double positionStride = 1e3;

    double safety = 0; //distance in nanometers we keep edeps from the edges of the sphere
    int nEdepsInSingleVolume = 1000;
    double inSphereStride = ((positionStride)-safety*2)/nEdepsInSingleVolume; //the insphere stride, is calcualted by subtracting the safety off each side. Then dividing by nEdepsInVolume
    
    for (int i = 0; i < 5; i++)
    {
        //Place a series of edeps along the z axis, each jumping 1 um 
        //Start at the edge of the box plus half a diameter length
        //Jump one diameter length every iteration
        for (double outputValue = startPosition; outputValue <= stopPosition; outputValue+=positionStride)
        {
            //Now we're in a sphere, place 40 edeps inside
            for (int j = 0; j < nEdepsInSingleVolume; j++)
            {
                //shifts us from center of sphere to edge, then increments by inSphereStride
                z = outputValue;//-(positionStride/2)+safety+(inSphereStride*j);
                //output
                pTrackOutputTree->Fill();
                indexEntry += 1;
            }
        }

        pEventIndexTree->Fill();
    }

    std::cout << indexEntry << std::endl;

    pTrackOutputFile->Write(0,TObject::kWriteDelete);
    pTrackOutputFile->Close();
}

void generateTestTracks::GenerateTestA6() // Generate a track where all the points are on the surface of the cube
{
    //Need to establish my variables
    double x,y,z,edep;
    long indexEntry = 0;

    //Set track information to known fixed values
    edep = double(2000)/double(3); //2000 eV to give 1 keV per micrometer in a 1 um diameter sphere

    //Create TFile
    auto pTrackOutputFile = new TFile("./testA6/testA6.root","RECREATE");

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
    
    double startPosition = -1.5e6+0.5e3;
    double stopPosition = 1.5e6-0.5e3;
    double positionStride = 1e3;

    /*
    NOTE: Remember that the current algorithm doesn't even check if particles are
    outside of the box in the Z-axis, because they're designed not to be. 

    That's why we're only checking 4 sides of the cube
    */

    //Side 1
    x = startPosition;
    for (z = startPosition; z <= stopPosition; z += positionStride)
    {
        for (y = startPosition; y <= stopPosition; y += positionStride)
        {
            //output
            pTrackOutputTree->Fill();
            indexEntry += 1;
        }
    }

    //Side 2
    x = stopPosition;
    for (z = startPosition; z <= stopPosition; z += positionStride)
    {
        for (y = startPosition; y <= stopPosition; y += positionStride)
        {
            //output
            pTrackOutputTree->Fill();
            indexEntry += 1;
        }
    }

    //Side 3
    y = startPosition;
    for (x = startPosition; x <= stopPosition; x += positionStride)
    {
        for (y = startPosition; y <= stopPosition; y += positionStride)
        {
            //output
            pTrackOutputTree->Fill();
            indexEntry += 1;
        }
    }

    //Side 4
    y = stopPosition;
    for (x = startPosition; x <= stopPosition; x += positionStride)
    {
        for (z = startPosition; z <= stopPosition; z += positionStride)
        {
            //output
            pTrackOutputTree->Fill();
            indexEntry += 1;
        }
    }
    
    pEventIndexTree->Fill();

    std::cout << indexEntry << std::endl;

    pTrackOutputFile->Write(0,TObject::kWriteDelete);
    pTrackOutputFile->Close();
}

void generateTestTracks::GenerateTestA7() // Similar to the track for A2, but more points in z
{
    //This is the track that we use for calculation of pi
    //To determine that our random number generator works well

    //Need to establish my variables
    double x,y,z,edep;
    long indexEntry = 0;

    //Set track information to known fixed values
    x = 0; y = 0; //X and Y along the axis of spheres nearest the center of the box
    edep = double(2000)/3; //2000 eV to give 1 keV per micrometer in a 1 um diameter sphere

    //Create TFile
    auto pTrackOutputFile = new TFile("./testA7/testA7.root","RECREATE");

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
    
    double startPosition = -0.1e6;
    double stopPosition = 0.1e6;
    double positionStride = 10;

   for (int i = 0; i < 1; i++)
    {
        //Place a series of edeps along the z axis, each jumping 1 um 
        //Start at the edge of the box plus half a diameter length
        //Jump one diameter length every iteration
        for (double outputValue = startPosition; outputValue < stopPosition; outputValue+=positionStride)
        {
            z = outputValue;
            pTrackOutputTree->Fill();
            indexEntry += 1;
        }

        pEventIndexTree->Fill();
    }

    std::cout << indexEntry << std::endl;

    pTrackOutputFile->Write(0,TObject::kWriteDelete);
    pTrackOutputFile->Close();
}

void generateTestTracks::GenerateTestA8() //Energy depositions at a point to check overflows
{
    //This test just shows whether numbers are overflowing at the 4-byte int level or not

    //To make this test even better. Copy this printf into VoxelConstrainedSphereMethodKernel::ScoreTrackInSphere
    //At the very end
    /*
    printf("X,Y,Z position: %f, %f, %f \n X,Y,Z index: %ld, %ld, %ld SphereHitIndex: %ld \n \n",inputTrack.x[trackIdInSphere[i]],inputTrack.y[trackIdInSphere[i]],inputTrack.z[trackIdInSphere[i]], xIndex,yIndex,zIndex,outputPair.volume[i]);
    */

    //Need to establish my variables
    double x,y,z,edep;
    long indexEntry = 0;

    //Set track information to known fixed values
    edep = double(20)/double(3); //to give 1 keV per micrometer in a 10 nm diameter sphere

    //Create TFile
    auto pTrackOutputFile = new TFile("./testA8/testA8.root","RECREATE");

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

    double boxedge = 1.5e6;

    //The overflow test works with 10 nm diameter spheres in a 3 mm side length box
    //Entering positions manually
    x = 5-boxedge; y = 5-boxedge; z = 5-boxedge; //nm for the (0,0,0) sphere 
    pTrackOutputTree->Fill();
    indexEntry += 1;

    x = 836465-boxedge; y = 71585-boxedge;  //one below the potential overflow
    pTrackOutputTree->Fill();
    indexEntry += 1;

    x = 836475-boxedge; y = 71585-boxedge;  //nm for the (836475,7158,0) sphere. Where it should overflow to negative numbers 
    pTrackOutputTree->Fill();
    indexEntry += 1;

    x = 836485-boxedge; y = 71585-boxedge;  //one above the potential overflow
    pTrackOutputTree->Fill();
    indexEntry += 1;

    x = 1672955-boxedge; y = 143166-boxedge; //where it would overflow with an unsigned 4 byte int 
    pTrackOutputTree->Fill();
    indexEntry += 1;

    x = 1672965-boxedge; y = 143166-boxedge; //one above where it would overflow with an unsigned 4 byte
    pTrackOutputTree->Fill();
    indexEntry += 1;

    x = 5-boxedge; y = 5-boxedge; z = 1.5e6-5;//this is just an astronomically large value. Test it with the print statement above to be sure it's working. 
    pTrackOutputTree->Fill();
    indexEntry += 1;

    pEventIndexTree->Fill();

    std::cout << indexEntry << std::endl;

    pTrackOutputFile->Write(0,TObject::kWriteDelete);
    pTrackOutputFile->Close();
}