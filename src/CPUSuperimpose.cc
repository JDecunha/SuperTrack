//SuperTrack
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

//SMatrix and SVector are the fastest
//ways to hold vectors and matrices in ROOT
typedef ROOT::Math::SVector<Double_t,3> SVector3;
#define VERBOSE 0
#define VERY_VERBOSE 1

using namespace std;

template <class datatype>
struct DataStruct
{

	datatype data; //an arbitrary object to either store a histogram or list of energy values
	std::unordered_map<string,string> information_strings; //a series of information strings

	DataStruct(datatype d) : data(d) {} //this is a constructor even though it doesn't look like it

	void save(TString output_folder)
	{
		/*
		So what gets saved?

		What do I want me filename to be?:
		[analyzed_filename]_[random_seed]_[nhistories]x[noversamples].root
		*/
		// if(typeid(datatype) == typeid(TH1F) {} is a way of checking the template name to change the save file

		//Generate the path and filename of the file to be saved

		string analyzed_filename = std::filesystem::path(string(information_strings["simulation_track_structure_filename"])).stem(); //extract just the filename of the tracks analyzed
		stringstream namestream;
		namestream << output_folder << "/" << analyzed_filename << "_" << information_strings["simulation_random_seed"] << "_" << information_strings["simulation_nhistories_analyzed"] << "x" << information_strings["simulation_noversamples"] << ".root";
		string output_filename = namestream.str();

		//Generate and open the .root File
		TFile savefile = TFile(output_filename.c_str(),"recreate");//recreate will create a new file, or overwrite if the file already exists

		//Put the objects to be saved in the .root into an array
		TObjArray savefile_array(0);
		//add the data
		//savefile_array.Add(&data);
		savefile.WriteObject(&data,"Data");

		//loop through the string and add to the object array
		std::unordered_map<string,string>::iterator it = information_strings.begin();
		while (it != information_strings.cend())
		{
			savefile_array.Add(new TNamed(TString(it->first),TString(it->second))); //this could be a memory leak, hmm
			it++;
		}
		savefile_array.Write();
		savefile.Close();
	}
};

TH1F score_lineal_voxel(TString filepath, float_t scoring_sphere_spacing, float_t scoring_sphere_diameter, Int_t nthreads, Int_t nSamples = 1, Long_t random_seed = time(NULL))
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

	float RAND_MAX_F = float(RAND_MAX);
	
	//the = sign captures everything in the enclosing function by value. Meaning it makes a process local copy.
	auto workItem = [=](std::tuple<Int_t,Int_t,Int_t,TString> input) 
	{
		//Open the file in each process and make a Tree Reader
		TFile f = TFile(get<3>(input));
		TTreeReader trackReader("Tracks", &f);
		trackReader.SetEntriesRange(get<0>(input),get<1>(input));
		TTreeReaderValue<double_t> x(trackReader, "x [nm]");
		TTreeReaderValue<double_t> y(trackReader, "y [nm]");
		TTreeReaderValue<double_t> z(trackReader, "z [nm]");
		TTreeReaderValue<double_t> edep(trackReader, "edep [eV]");

		cout << "thread #: " << get<2>(input) << " starting at: " << to_string(get<0>(input)) << endl;

		//Get the voxel side length
		TNamed *voxelSideLengthName;
		double_t voxelSideLength;
		f.GetObject("Voxel side length [mm]",voxelSideLengthName);
		voxelSideLength = 1e6*atof(voxelSideLengthName->GetTitle()); //Get voxel side length and convert to nm
		double_t scoring_square_half_length = voxelSideLength/2;

		//Initialize the geometry
		long long int index;
		int num_spheres_linear = TMath::Ceil(((scoring_square_half_length*2)/scoring_sphere_spacing)); //this is how many spheres there will be in a line
		if  (num_spheres_linear % 2 == 0) 
		{
			#if VERBOSE == 2
			cout << "even, adding 1 extra sphere" << endl; 
			#endif
			num_spheres_linear += 1; //This is so there is always a central, 0,0,0 sphere
		} 
		long long int num_spheres_total = TMath::Power((num_spheres_linear),3);
		float_t top_sphere_offset = -(((float(num_spheres_linear))/2)-0.5)*scoring_sphere_spacing;//So the furthest sphere away in X,Y, or Z will be number of spheres plus the half center sphere away from the center

		#if VERBOSE == 2
		cout << "Number of spheres in a line: " << num_spheres_linear << endl;
		cout << "Total number of spheres in a cube: " << num_spheres_total << endl;
		cout << "largest sphere offset from central sphere: " << top_sphere_offset << endl;
		#endif

		//Initialize the histogram
		TH1F lineal_histogram = TH1F("Lineal energy histogram", "y*f(y)", 200, -2,1);
		BinLogX(&lineal_histogram); //transform the bins to logarithmic
		std::unordered_map<int,double> energy_map; //define an unordered dictionary to hold edep values and associated volume

		//Initialize the random number generator. Append the thread ID to the current time
		Long_t random_seed = std::stol(std::to_string(random_seed) + std::to_string(get<2>(input)));
		srand(random_seed); //TODO: test that my threads are getting different random numbers (I know they're getting different seeds, but does this work?)

		//Initalize the vectors and matrices we're going to use in the looping
		SVector3 particle_position, position_shifts, position_difference_from_nearest_sphere,position_indices;


		//Currently we haven't implemented supersampling. Have to think about how we recombine histograms and what not
		//Start a loop for each time you sample the track
		for (int i = 0; i < nSamples; i++)
		{

			//set the position shifts. No shift in Z axis, because we only want to shift on the X-Y surface of the box
			position_shifts[0] = ((rand())*scoring_square_half_length*2/(RAND_MAX))-scoring_square_half_length;
			position_shifts[1] = ((rand())*scoring_square_half_length*2/(RAND_MAX))-scoring_square_half_length;
			position_shifts[2] = 0;

			while (trackReader.Next())
			{
				//Set the particle's position into a vector
				particle_position[0] = *x+position_shifts[0];
				particle_position[1] = *y+position_shifts[1];
				particle_position[2] = *z;

				//Check if inside box
					if (abs(particle_position[0]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2) && abs(particle_position[1]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2) && abs(particle_position[2]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2)) // if inside box
					{

						//Convert from x,y,z to index position
						position_indices[0] = TMath::Nint((particle_position[0]-top_sphere_offset)/scoring_sphere_spacing);
						position_indices[1] = TMath::Nint((particle_position[1]-top_sphere_offset)/scoring_sphere_spacing);
						position_indices[2] = TMath::Nint((particle_position[2]-top_sphere_offset)/scoring_sphere_spacing);

						//Figure out if x,y,z coordinate is within sphere
						//top_sphere_offset+(x_index*scoring_sphere_spacing) should give you the center of the sphere closest to your coordinate
						position_difference_from_nearest_sphere = position_indices*scoring_sphere_spacing;
						position_difference_from_nearest_sphere = position_difference_from_nearest_sphere+top_sphere_offset-particle_position;

						if(ROOT::Math::Mag(position_difference_from_nearest_sphere) <= (scoring_sphere_diameter/2))
						{
							//Okay you are inside the sphere
							index = position_indices[0] + (position_indices[1]*(num_spheres_linear)) + position_indices[2]*TMath::Power((num_spheres_linear),2); //Keep in mind that for the index it starts counting at zero
							double_t lineal_energy = *edep/((2./3.)*scoring_sphere_diameter); //this should be ev/nm which is same a kev/um
							energy_map[index] += lineal_energy; 
						} 
				 }
			
		  	}
		//Has finished iterating over current event. Output data
		trackReader.Restart(); //set iterator back to the beginning
		std::unordered_map<int,double_t>::iterator it = energy_map.begin();
		while (it != energy_map.cend())
		{
			lineal_histogram.Fill(it->second);
			it++;
		}
		energy_map.clear(); //empty the map for next event

		}

	  return lineal_histogram;

	};


	// Create the pool of workers
  ROOT::TProcessExecutor workers(nthreads);
  //Process the jobs and get a vector of the output
  std::vector<TH1F> process_output = workers.Map(workItem, perthread_input_arguments);

   //THIS IS SO JANKY
   //But according to the CERN forms this is the 'best' way 
   //1.) make a copy of the TH1F and point to it
   //2.) dynamically allocate a list
   //3.) use the new TH1F to merge the list
   //4.) Take the TH1F pointer and store the value it points to in stack memory
   //5.) Delete the list so I don't leak memory
   //6.) Return the local merged TH1F
   TH1F* h = (TH1F*)process_output[0].Clone();
   TList *list = new TList;
   for (TObject & hist : process_output)
   {
   	list->Add(&hist);
   }
   h->Merge(list);
   TH1F h2 = *h;
   delete list; //Delete the list so I don't leak memory

 	 //PMF_to_PDF(&h2);
	 //Prepare_for_Semilog(&h2);

   return h2;

}