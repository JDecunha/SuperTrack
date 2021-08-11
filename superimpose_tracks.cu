#include "TFile.h"
#include "THStack.h"
#include "TPad.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TROOT.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TTree.h"
#include "TMath.h"
#include "TEntryList.h"
#include "ROOT/TProcessExecutor.hxx"
#include "Math/SMatrix.h"
#include "Math/LorentzVector.h"
#include "Math/Vector4Dfwd.h"
#include <unordered_map>
#include <vector>
#include <iterator>
#include <tuple>
#include <filesystem>
#include "include/testCUDA.cuh"
#include <stdio.h>


//SMatrix and SVector are the fastest
//ways to hold vectors and matrices in ROOT
typedef ROOT::Math::SMatrix<Double_t,3> SMatrix33;
typedef ROOT::Math::SVector<Double_t,3> SVector3;
typedef ROOT::Math::SVector<int,3> SIntegerVector3;
#define VERBOSE 0
#define VERY_VERBOSE 1

using namespace std;

//GLOBAL VARIABLES
//TTree *microdosimetry;

//USEFUL COMMANDS
//from command line: Tree->Show(entry number) is useful
//also Tree->Print() shows all the info but no specific entries
// .x file.c+ executes and compiles in root

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

int determine_number_of_events(TTree *tree)
{
	//TODO: see if I can either 1.) get the previous branch address and reset it after I set it here 
	//2.) OR just find a way to access value without having to set the branch address at all. Becuase this seems like a poor way of accessing a variable.
	int number_of_edeps = tree->GetEntries();
	int number_of_events;
	tree->SetBranchAddress("eventID",&number_of_events);
	tree->GetEntry(number_of_edeps-1);
	tree->ResetBranchAddresses();

	return number_of_events+1;
}

void transform_microdosimetry_file_to_4_vector(TString filepath, bool process_all, Int_t nevents=0)
{

	TFile f(filepath);
	if (process_all == true) 
	{
		TTree* microdosimetry;
		f.GetObject("microdosimetry",microdosimetry);
		nevents = determine_number_of_events(microdosimetry);
	}

	string analyzed_filename = std::filesystem::path(string(filepath)).stem(); //extract just the filename of the tracks analyzed
	string path_to_io_folder = std::filesystem::path(string(filepath)).parent_path(); //extract the path of the tracks analyzed
	
	stringstream namestream;
	namestream << path_to_io_folder << "/" << analyzed_filename << "_4vectortransformed_" << to_string(nevents) << ".root";
	string output_filename = namestream.str();

	//Generate and open the .root File
	TFile savefile = TFile(output_filename.c_str(),"recreate");//recreate will create a new file, or overwrite if the file already exists
	
	
	//okay, time to read the original Tree
	TTreeReader microdosimetryreader("microdosimetry", &f);
	TTreeReaderValue<Double_t> x(microdosimetryreader, "x");
	TTreeReaderValue<Double_t> y(microdosimetryreader, "y");
	TTreeReaderValue<Double_t> z(microdosimetryreader, "z");
	TTreeReaderValue<Int_t> eventID(microdosimetryreader, "eventID");
	TTreeReaderValue<Double_t> edep(microdosimetryreader, "totalEnergyDeposit");

	//Make a new TTree
	//with a 4 vector TBranch
	//and an eventID TBranch
	TTree modifiedtree("microdosimetry","microdosimetry");
	std::vector<ROOT::Math::XYZTVector> fourvector;
	Int_t current_event_number = 0;
	auto fourvectorbranch = modifiedtree.Branch("XYZEdepVector",&fourvector,32000,0); //We want 0 splitlevel because we always access X,Y,Z,Edep together, so write them contigiously in memory
	auto eventIDbranch = modifiedtree.Branch("eventID",&current_event_number,32000,0);


	while (microdosimetryreader.Next())
	{
		if (current_event_number == *eventID)
		{
			fourvector.push_back(ROOT::Math::XYZTVector(*x,*y,*z,*edep));
		}
		else
		{
			modifiedtree.Fill();
			fourvector.clear();
			if(*eventID >= nevents-1)
			{
				break;
			}
			else
			{
			fourvector.push_back(ROOT::Math::XYZTVector(*x,*y,*z,*edep)); //write the next element to a fresh fourvector
			current_event_number = *eventID;
			}
		}
	}
	if (process_all == true){modifiedtree.Fill();} //if you're processing all, the loop stops before the last fill. So call it here.

	f.Close();
	savefile.Write("", TObject::kOverwrite);//Object(&modifiedtree,"microdosimetry");
	savefile.Close();
	
	
}

void BinLogX(TH1* h)
{
	///Very useful function from the ROOT forums
	///If you enter your axes in Log10 notation (i.e. -3,0) in the initial TH1 constructor
	///This will reformat your bins logarithmically for you

   TAxis *axis = h->GetXaxis();
   int bins = axis->GetNbins();

   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];

   for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);

   }
   axis->Set(bins, new_bins);
   delete[] new_bins;
} 
void BinLogXMultithread(std::shared_ptr<TH1F> h)
{
	///Very useful function from the ROOT forums
	///If you enter your axes in Log10 notation (i.e. -3,0) in the initial TH1 constructor
	///This will reformat your bins logarithmically for you

   TAxis *axis = h->GetXaxis();
   int bins = axis->GetNbins();

   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];

   for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);

   }
   axis->Set(bins, new_bins);
   delete[] new_bins;
} 

void PMF_to_PDF(TH1* h)
{
	//This converts a probability mass function to a probability density function
	int length = h->GetNbinsX();
	float_t normalization = 0;
	for (int i = 0; i < length; i++) //bins per bin width
	{
		auto value = h->GetBinContent(i);
		auto low_edge = h->GetBinLowEdge(i);
		auto high_edge = h->GetBinLowEdge(i+1);
		auto per_width = value/(high_edge-low_edge);
		h->SetBinContent(i,per_width);

		normalization += value;
	}
	for (int i = 0; i < length; i++) //normalize
	{
		auto value = h->GetBinContent(i);
		auto normalized = value/normalization;
		h->SetBinContent(i,normalized);
	}
}

void Prepare_for_Semilog(TH1* h)
{
	//This takes a PDF and normalizes it by an extra factor of y (i.e. to give y*f(y))to preserve the graphical properties
	//of a PDF on a semilog axis
	int length = h->GetNbinsX();
	for (int i = 0; i < length; i++) //bins per bin width
	{
		auto value = h->GetBinContent(i);
		auto low_edge = h->GetBinLowEdge(i);
		auto high_edge = h->GetBinLowEdge(i+1);
		auto mid_value = (high_edge+low_edge)/2;
		h->SetBinContent(i,value*mid_value);
	}
}

SMatrix33 uniform_random_rotation_matrix(float x0,float x1,float x2)
{
	//input: x0,x1, and x2 are random numbers sampled between 0 and 1
	//output: a uniformly sampled random 3D rotation matrix

	//
	// Though I trust the book this is derived from I should do some benchmarking just to verify
	// that my implementation is working correctly.
	// Benchmarking would basically involve generating histograms of x,y,z locations on the sphere and checking they line up
	//
	// This version of the random rotation looks to be working
	// at least the determinant and the transpose of the matrix are the same
	// as expected for a rotation matrix
	//
	
	//Adapted from Graphics Gems III Fast Random Rotation Matrices, Section by James Avro, appendix includes C code algorithm
	Double_t theta,phi,z;
	Double_t r,st,ct;
	Double_t Sx,Sy;
	Double_t Vx, Vy, Vz;

	//Normalize the input random numbers to new values
	theta = x0 * TMath::Pi() * 2;
	phi = x1 * TMath::Pi() * 2;
	z = x2 * 2;

	r = TMath::Sqrt(z);
	Vx = TMath::Sin(phi)*r;
	Vy = TMath::Cos(phi)*r;
	Vz = TMath::Sqrt(2-z);

	st = TMath::Sin(theta);
	ct = TMath::Cos(theta);
	Sx = Vx * ct - Vy * st;
	Sy = Vx * st + Vy * ct;

	Double_t flat_matrix[9] = {Vx*Sx-ct,Vx*Sy-st,Vx*Vz,Vy*Sx+st,Vy*Sy-ct,Vy*Vz,Vz*Sx,Vz*Sy,1-z};

	return SMatrix33(flat_matrix,9);
}

void uniform_random_rotation_matrix_optimized(float x0,float x1,float x2, SMatrix33* matrix)
{
	//input: x0,x1, and x2 are random numbers sampled between 0 and 1
	//output: a uniformly sampled random 3D rotation matrix

	//
	// Though I trust the book this is derived from I should do some benchmarking just to verify
	// that my implementation is working correctly.
	// Benchmarking would basically involve generating histograms of x,y,z locations on the sphere and checking they line up
	//
	// This version of the random rotation looks to be working
	// at least the determinant and the transpose of the matrix are the same
	// as expected for a rotation matrix
	//
	
	//cout << " x: " << x0 << " y: " << x1 << " z: " << x2 << endl;

	//Adapted from Graphics Gems III Fast Random Rotation Matrices, Section by James Avro, appendix includes C code algorithm
	Double_t theta,phi,z;
	Double_t r,st,ct;
	Double_t Sx,Sy;
	Double_t Vx, Vy, Vz;

	//Normalize the input random numbers to new values
	theta = x0 * TMath::Pi() * 2;
	phi = x1 * TMath::Pi() * 2;
	z = x2 * 2;

	r = TMath::Sqrt(z);
	Vx = TMath::Sin(phi)*r;
	Vy = TMath::Cos(phi)*r;
	Vz = TMath::Sqrt(2-z);

	st = TMath::Sin(theta);
	ct = TMath::Cos(theta);
	Sx = Vx * ct - Vy * st;
	Sy = Vx * st + Vy * ct;

	double matrixvals[9] = {Vx*Sx-ct,Vx*Sy-st,Vx*Vz,Vy*Sx+st,Vy*Sy-ct,Vy*Vz,Vz*Sx,Vz*Sy,1-z};

	matrix->SetElements(matrixvals,matrixvals+9);
}

TEntryList* create_entry_list(TFile* f, Int_t nevents)
{

	//TFile *f = TFile::Open(filepath);
	TTreeReader microdosimetryreader("microdosimetry", f);
	TTreeReaderValue<Int_t> eventID(microdosimetryreader, "eventID");
	TEntryList* microdosimetry_entry_list = new TEntryList();

	while (microdosimetryreader.Next())
	{
		if (*eventID < nevents)
		{
			microdosimetry_entry_list->Enter(microdosimetryreader.GetCurrentEntry());
		}
		else {break;}
	}

	return microdosimetry_entry_list;

}

/*void score_lineal_histogram_multithreaded(TString filepath, float_t scoring_square_half_length, float_t scoring_sphere_spacing, float_t scoring_sphere_diameter, int track_oversample_factor, float_t CPE_range, Long_t random_seed = time(NULL),Int_t nhistoriestoanalyze = 0)
{

	
	DESCRIPTION OF THE ALGORITHM:
	1.) Rotation the edep point by a random rotation matrix.
	2.) Shift the edep point by a range
	3.) Oversample tracks by track_oversample_factor number of times
	4.) Do not try to place the track in a volume if the track is outside of the box region
	5.) Calculate lineal energy from that
	
	int nthreads = 4;
	ROOT::EnableImplicitMT(nthreads);

	//Open the file, retrieve the tree, and determine the number of histories
	TFile *f = TFile::Open(filepath);
	TTree *microdosimetry;
	f->GetObject("microdosimetry",microdosimetry);

	//TODO: see if setting microdosimetry's entry list outside of TTreeProcessor works as well
	auto entry_list = create_entry_list(f,100000);
	ROOT::TTreeProcessorMT tp(*microdosimetry,*entry_list);

	//TODO: See if I can modify this object's bins here and not in the threads
	//If I can't modify the histogram here
	//I could create a regular TH1F, change its bins, and then assign the TThreaded object with a copy constructor
  	ROOT::TThreadedObject<TH1F> histogram("histogram","y*f(y)",  200, -2, 1.5);

  	//Initialize the geometry
	long long int index;
	int num_spheres_linear = TMath::Ceil(((scoring_square_half_length*2)/scoring_sphere_spacing)); //this is how many spheres there will be in a line
	if  (num_spheres_linear % 2 == 0) {cout << "even, adding 1 extra sphere" << endl; num_spheres_linear += 1;} //This is so there is always a central, 0,0,0 sphere
	long long int num_spheres_total = TMath::Power((num_spheres_linear),3);
	float_t top_sphere_offset = -(((float(num_spheres_linear))/2)-0.5)*scoring_sphere_spacing;//So the furthest sphere away in X,Y, or Z will be number of spheres plus the half center sphere away from the center

	// Define the function that will process a subrange of the tree.
	// The function must receive only one parameter, a TTreeReader,
	// and it must be thread safe. To enforce the latter requirement,
	// TThreadedObject histograms will be used.
	auto myFunction = [&](TTreeReader &myReader) 
	{
	 TTreeReaderValue<Double_t> xRV(myReader, "x");
	 TTreeReaderValue<Double_t> yRV(myReader, "y");
	 TTreeReaderValue<Double_t> zRV(myReader, "z");
	 //Try to make this and everything else thread safe too
	 SVector3 particle_position = SVector3();
	 //TODO: make this unordered map thread safe, and move it outside the function
	 std::unordered_map<int,double> energy_map; //define an unordered dictionary to hold edep values and associated volume
	 TTreeReaderValue<Int_t> eventIDRV(myReader, "eventID");
	 cout << "Line invoked" << endl;
	 cout << std::this_thread::get_id() << endl;

	 // For performance reasons, a copy of the pointer associated to this thread on the
	 // stack is used
	 auto thread_local_hist = histogram.Get();
	 BinLogXMultithread(thread_local_hist); //transform the bins to logarithmic

	 //Set values to initialize loop
	 int k = 0;
	 int inside = 0; //a counter for the number of edeps within the box AND within a sphere
	 int outside = 0; //a counter for the number of edeps within the box BUT outside of a sphere
	 srand(random_seed); //Give current time as a seed to random number generator

	 while (myReader.Next()) 
	 {
	    particle_position[0] = *xRV;
	    particle_position[1] = *yRV;
	    particle_position[2] = *zRV;

	    auto eventID = *eventIDRV;

		#if VERBOSE == 1
		//cout << particle_position[0] << endl;
		#endif
	 }
	    
	};

	// Launch the parallel processing of the tree
	tp.Process(myFunction);

	// Use the TThreadedObject::Merge method to merge the thread private histograms
	// into the final result
	//auto histogram_merged = histogram.Merge();

	return 0;
}*/

TH1F score_lineal_histogram_multithreaded_explicit(TString filepath, float_t scoring_square_half_length, float_t scoring_sphere_spacing, float_t scoring_sphere_diameter, float_t CPE_range,Int_t nthreads, Long_t random_seed = time(NULL),Long64_t nhistoriestoanalyze = 0)
{
	//open the file, retrive the tree
	TFile f = TFile(filepath);
	TTree *microdosimetry;
	f.GetObject("microdosimetry",microdosimetry);

	//Now, loop over the Tree cluster by cluster,
	//and determine the entry and exit points for each TreeReader
	Long64_t nentries = microdosimetry->GetEntries();
	if (nhistoriestoanalyze < nentries && nhistoriestoanalyze != 0) {nentries = nhistoriestoanalyze;} //only if the number of histories you have requested is smaller than the file, then analyze that many. If 0, then analyze all.
	auto clusteriterator = microdosimetry->GetClusterIterator(0);
	std::vector<std::tuple<Int_t,Int_t,Int_t,TString>> input_arguments_multithreading;
	Long64_t start_entry_val = 0;
	Long64_t end_entry_val = 0;

	for (Int_t i = 1; i <= nthreads; i++)
	{
		while((end_entry_val = clusteriterator()) < (nentries*i)/nthreads) {}
		input_arguments_multithreading.push_back(std::make_tuple(start_entry_val,end_entry_val,i,filepath));
		//cout << "thread: " << i << " start: " << start_entry_val << " end: " << end_entry_val << endl;
		start_entry_val = end_entry_val + 1;
	}

	
	f.Close(); //I can close this file I'm done with it

	Long_t random_seed_base = time(NULL);

	float RAND_MAX_F = float(RAND_MAX);
	
	//the = sign captures everything in the enclosing function by value. Meaning it makes a process local copy.
	auto workItem = [=](std::tuple<Int_t,Int_t,Int_t,TString> input) 
	{
		//Open the file in each process and make a Tree Reader
		TFile f = TFile(get<3>(input));
		TTreeReader microdosimetryreader("microdosimetry", &f);
		microdosimetryreader.SetEntriesRange(get<0>(input),get<1>(input));
		TTreeReaderArray<ROOT::Math::XYZTVector> XYZEdep(microdosimetryreader, "XYZEdepVector");
		TTreeReaderValue<Int_t> eventID(microdosimetryreader, "eventID");

		cout << "thread #: " << get<2>(input) << " starting at: " << to_string(get<0>(input)) << endl;

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
		std::unordered_map<int,double> energy_map; //define an unordered dictionary to hold edep values and associated volume

		#if VERBOSE == 2
		cout << "Number of spheres in a line: " << num_spheres_linear << endl;
		cout << "Total number of spheres in a cube: " << num_spheres_total << endl;
		cout << "largest sphere offset from central sphere: " << top_sphere_offset << endl;
		cout << "Total number of entries in file: " << microdosimetry->GetEntries() << endl;
		#endif

		//Initialize the histogram
		TH1F lineal_histogram = TH1F("Lineal energy histogram", "y*f(y)", 200, -2,1);
		BinLogX(&lineal_histogram); //transform the bins to logarithmic

		//Initialize the random number generator. Append the thread ID to the current time
		Long_t random_seed = std::stol(std::to_string(random_seed_base) + std::to_string(get<2>(input)));
		//TODO: test that my threads are getting different random numbers (I know they're getting different seeds, but does this work?)
		srand(random_seed);

		//Initalize the vectors and matrices we're going to use in the looping
		SVector3 particle_position, position_shifts, position_rotated, position_shifted, position_difference_from_nearest_sphere,position_indices;
		//SIntegerVector3 position_indices;
		SMatrix33 rotation_matrix;

		while (microdosimetryreader.Next())
		{
			//set the position shifts
			position_shifts[0] = ((rand())*CPE_range*2/(RAND_MAX))-CPE_range;
			position_shifts[1] = ((rand())*CPE_range*2/(RAND_MAX))-CPE_range;
			position_shifts[2] = ((rand())*CPE_range*2/(RAND_MAX))-CPE_range;

			//refresh the rotation matrix
			uniform_random_rotation_matrix_optimized(rand()/RAND_MAX_F,rand()/RAND_MAX_F,rand()/RAND_MAX_F,&rotation_matrix);

			for (const ROOT::Math::XYZTVector & fourvector : XYZEdep)//(iterator = *XYZEdep->begin(); iterator != *XYZEdep->end(); iterator++)
			{
				//TODO: Compare these two methods
				//Because the ones below might not invoke the constructor
				// where this might
				//particle_position = fourvector.Vect();
				
				//My feeling is that this is slower than it needs to be because it goes component by component
				//See if I'm able to fill the 3 vector with an iterator or something a little faster.
				particle_position[0] = fourvector.X();
				particle_position[1] = fourvector.Y();
				particle_position[2] = fourvector.Z();

				if(fourvector.T() > 0) //check that there has been an energy deposition
				{
					//Transform particle posityion by the random rotation. Have to do this before shifting so the origin remains the same
					position_rotated = rotation_matrix*particle_position;

					//Shift the rotated particle position by a random x,y,z shift
					position_shifted = position_rotated+position_shifts;

					//Check if inside box
					if (abs(position_shifted[0]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2) && abs(position_shifted[1]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2) && abs(position_shifted[2]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2)) // if inside box
					{

						//Convert from x,y,z to index position
						position_indices[0] = TMath::Nint((position_shifted[0]-top_sphere_offset)/scoring_sphere_spacing);
						position_indices[1] = TMath::Nint((position_shifted[1]-top_sphere_offset)/scoring_sphere_spacing);
						position_indices[2] = TMath::Nint((position_shifted[2]-top_sphere_offset)/scoring_sphere_spacing);

						//Figure out if x,y,z coordinate is within sphere
						//top_sphere_offset+(x_index*scoring_sphere_spacing) should give you the center of the sphere closest to your coordinate
						position_difference_from_nearest_sphere = position_indices*scoring_sphere_spacing;
						position_difference_from_nearest_sphere = position_difference_from_nearest_sphere+top_sphere_offset-position_shifted;

						if(ROOT::Math::Mag(position_difference_from_nearest_sphere) <= (scoring_sphere_diameter/2))
						{
							//Okay you are inside the sphere
							index = position_indices[0] + (position_indices[1]*(num_spheres_linear)) + position_indices[2]*TMath::Power((num_spheres_linear),2); //Keep in mind that for the index it starts counting at zero
							double_t lineal_energy = fourvector.T()/((2./3.)*scoring_sphere_diameter); //this should be ev/nm which is same a kev/um
							energy_map[index] += lineal_energy; 
							//energy_map[index] += fourvector.T();
							//inside += 1;

						} //else {outside += 1;} 
					}
				}

			}
			//Has finished iterating over current event. Output data
			std::unordered_map<int,double_t>::iterator it = energy_map.begin();
 			while (it != energy_map.cend())
 			{
 				lineal_histogram.Fill(it->second);
 				it++;
 			}
 			energy_map.clear(); //empty the map for next event
		}

		PMF_to_PDF(&lineal_histogram);
		Prepare_for_Semilog(&lineal_histogram);

  	return lineal_histogram;
	};

   // Create the pool of workers
   ROOT::TProcessExecutor workers(nthreads);
   //Process the jobs and get a vector of the output
   std::vector<TH1F> process_output = workers.Map(workItem, input_arguments_multithreading);

   //THIS IS SO JANKY
   //I should consider posting on the CERN Root forums to get a better solution
   //in steps
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

   return h2;
}

TH1F debug_code_verify_rotation_shifts(TString filepath, float_t scoring_square_half_length, float_t scoring_sphere_spacing, float_t scoring_sphere_diameter, float_t CPE_range,Int_t nthreads, Long_t random_seed = time(NULL),Long64_t nhistoriestoanalyze = 0, int axis = 0)
{
	//open the file, retrive the tree
	TFile f = TFile(filepath);
	TTree *microdosimetry;
	f.GetObject("microdosimetry",microdosimetry);

	//Now, loop over the Tree cluster by cluster,
	//and determine the entry and exit points for each TreeReader
	Long64_t nentries = microdosimetry->GetEntries();
	if (nhistoriestoanalyze < nentries) {nentries = nhistoriestoanalyze;} //only if the number of histories you have requested is smaller than the file, then analyze that many. (i.e. you can't analyze histories you don't have)
	auto clusteriterator = microdosimetry->GetClusterIterator(0);
	std::vector<std::tuple<Int_t,Int_t,Int_t,TString>> input_arguments_multithreading;
	Long64_t start_entry_val = 0;
	Long64_t end_entry_val = 0;

	for (Int_t i = 1; i <= nthreads; i++)
	{
		while((end_entry_val = clusteriterator()) < (nentries*i)/nthreads) {}
		input_arguments_multithreading.push_back(std::make_tuple(start_entry_val,end_entry_val,i,filepath));
		//cout << "thread: " << i << " start: " << start_entry_val << " end: " << end_entry_val << endl;
		start_entry_val = end_entry_val + 1;
	}

	
	f.Close(); //I can close this file I'm done with it

	Long_t random_seed_base = time(NULL);
	
	//the = sign captures everything in the enclosing function by value. Meaning it makes a process local copy.
	auto workItem = [=](std::tuple<Int_t,Int_t,Int_t,TString> input) 
	{
		//Open the file in each process and make a Tree Reader
		TFile f = TFile(get<3>(input));
		TTreeReader microdosimetryreader("microdosimetry", &f);
		microdosimetryreader.SetEntriesRange(get<0>(input),get<1>(input));
		TTreeReaderArray<ROOT::Math::XYZTVector> XYZEdep(microdosimetryreader, "XYZEdepVector");
		TTreeReaderValue<Int_t> eventID(microdosimetryreader, "eventID");

		cout << "thread #: " << get<2>(input) << " starting at: " << to_string(get<0>(input)) << endl;

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
		std::unordered_map<int,double> energy_map; //define an unordered dictionary to hold edep values and associated volume

		#if VERBOSE == 2
		cout << "Number of spheres in a line: " << num_spheres_linear << endl;
		cout << "Total number of spheres in a cube: " << num_spheres_total << endl;
		cout << "largest sphere offset from central sphere: " << top_sphere_offset << endl;
		cout << "Total number of entries in file: " << microdosimetry->GetEntries() << endl;
		#endif

		//Initialize the histogram
		TH1F lineal_histogram = TH1F("Lineal energy histogram", "y*f(y)", 200, -1.1,1.1);

		//Initialize the random number generator. Append the thread ID to the current time
		Long_t random_seed = std::stol(std::to_string(random_seed_base) + std::to_string(get<2>(input)));
		//TODO: test that my threads are getting different random numbers (I know they're getting different seeds, but does this work?)
		srand(random_seed);

		//Initalize the vectors and matrices we're going to use in the looping
		SVector3 particle_position, position_shifts, position_rotated, position_shifted, position_difference_from_nearest_sphere,position_indices;
		//SIntegerVector3 position_indices;
		SMatrix33 rotation_matrix;

		float RAND_MAX_F = float(RAND_MAX);

		while (microdosimetryreader.Next())
		{
			//refresh the rotation matrix
			uniform_random_rotation_matrix_optimized(rand()/RAND_MAX_F,rand()/RAND_MAX_F,rand()/RAND_MAX_F,&rotation_matrix);
			//cout << rotation_matrix(0,0) << " " << rotation_matrix(0,1) << " " << rotation_matrix(0,2) << endl; 
			//cout << rotation_matrix(1,0) << " " << rotation_matrix(1,1) << " " << rotation_matrix(1,2) << endl; 
			//cout << rotation_matrix(2,0) << " " << rotation_matrix(2,1) << " " << rotation_matrix(2,2) << endl; 

			for (const ROOT::Math::XYZTVector & fourvector : XYZEdep)//(iterator = *XYZEdep->begin(); iterator != *XYZEdep->end(); iterator++)
			{
				//TODO: Compare these two methods
				//Because the ones below might not invoke the constructor
				// where this might
				//particle_position = fourvector.Vect();
				
				//My feeling is that this is slower than it needs to be because it goes component by component
				//See if I'm able to fill the 3 vector with an iterator or something a little faster.
				particle_position[0] = 0;
				particle_position[1] = 0;
				particle_position[2] = 1;

				if(fourvector.T() > 0) //check that there has been an energy deposition
				{
					//Transform particle posityion by the random rotation. Have to do this before shifting so the origin remains the same
					position_rotated = rotation_matrix*particle_position;
					//float_t x_rotated = rotation_matrix(0,0)*particle_position[0]+rotation_matrix(0,1)*particle_position[1]+rotation_matrix(0,2)*particle_position[2];
					lineal_histogram.Fill(position_rotated[axis]);
				}

			}
			
		}
  	return lineal_histogram;
	};

   // Create the pool of workers
   ROOT::TProcessExecutor workers(nthreads);
   //Process the jobs and get a vector of the output
   std::vector<TH1F> process_output = workers.Map(workItem, input_arguments_multithreading);

   //THIS IS SO JANKY
   //I should consider posting on the CERN Root forums to get a better solution
   //in steps
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

   return h2;
}

TH1F debug_code_verify_position_shifts(TString filepath, float_t scoring_square_half_length, float_t scoring_sphere_spacing, float_t scoring_sphere_diameter, float_t CPE_range,Int_t nthreads, Long_t random_seed = time(NULL),Long64_t nhistoriestoanalyze = 0,int axis = 0)
{
	//open the file, retrive the tree
	TFile f = TFile(filepath);
	TTree *microdosimetry;
	f.GetObject("microdosimetry",microdosimetry);

	//Now, loop over the Tree cluster by cluster,
	//and determine the entry and exit points for each TreeReader
	Long64_t nentries = microdosimetry->GetEntries();
	if (nhistoriestoanalyze < nentries) {nentries = nhistoriestoanalyze;} //only if the number of histories you have requested is smaller than the file, then analyze that many. (i.e. you can't analyze histories you don't have)
	auto clusteriterator = microdosimetry->GetClusterIterator(0);
	std::vector<std::tuple<Int_t,Int_t,Int_t,TString>> input_arguments_multithreading;
	Long64_t start_entry_val = 0;
	Long64_t end_entry_val = 0;

	for (Int_t i = 1; i <= nthreads; i++)
	{
		while((end_entry_val = clusteriterator()) < (nentries*i)/nthreads) {}
		input_arguments_multithreading.push_back(std::make_tuple(start_entry_val,end_entry_val,i,filepath));
		//cout << "thread: " << i << " start: " << start_entry_val << " end: " << end_entry_val << endl;
		start_entry_val = end_entry_val + 1;
	}

	
	f.Close(); //I can close this file I'm done with it

	Long_t random_seed_base = time(NULL);
	
	//the = sign captures everything in the enclosing function by value. Meaning it makes a process local copy.
	auto workItem = [=](std::tuple<Int_t,Int_t,Int_t,TString> input) 
	{
		//Open the file in each process and make a Tree Reader
		TFile f = TFile(get<3>(input));
		TTreeReader microdosimetryreader("microdosimetry", &f);
		microdosimetryreader.SetEntriesRange(get<0>(input),get<1>(input));
		TTreeReaderArray<ROOT::Math::XYZTVector> XYZEdep(microdosimetryreader, "XYZEdepVector");
		TTreeReaderValue<Int_t> eventID(microdosimetryreader, "eventID");

		cout << "thread #: " << get<2>(input) << " starting at: " << to_string(get<0>(input)) << endl;

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
		std::unordered_map<int,double> energy_map; //define an unordered dictionary to hold edep values and associated volume

		#if VERBOSE == 2
		cout << "Number of spheres in a line: " << num_spheres_linear << endl;
		cout << "Total number of spheres in a cube: " << num_spheres_total << endl;
		cout << "largest sphere offset from central sphere: " << top_sphere_offset << endl;
		cout << "Total number of entries in file: " << microdosimetry->GetEntries() << endl;
		#endif

		//Initialize the histogram
		TH1F lineal_histogram = TH1F("Lineal energy histogram", "y*f(y)", 200, -6e6,6e6);

		//Initialize the random number generator. Append the thread ID to the current time
		Long_t random_seed = std::stol(std::to_string(random_seed_base) + std::to_string(get<2>(input)));
		//TODO: test that my threads are getting different random numbers (I know they're getting different seeds, but does this work?)
		srand(random_seed);

		//Initalize the vectors and matrices we're going to use in the looping
		SVector3 particle_position, position_shifts, position_rotated, position_shifted, position_difference_from_nearest_sphere,position_indices;
		//SIntegerVector3 position_indices;
		SMatrix33 rotation_matrix;

		while (microdosimetryreader.Next())
		{
			//set the position shifts
			position_shifts[0] = ((rand())*CPE_range*2/(RAND_MAX))-CPE_range;
			position_shifts[1] = ((rand())*CPE_range*2/(RAND_MAX))-CPE_range;
			position_shifts[2] = ((rand())*CPE_range*2/(RAND_MAX))-CPE_range;

			lineal_histogram.Fill(position_shifts[axis]);
		}


  	return lineal_histogram;
	};

   // Create the pool of workers
   ROOT::TProcessExecutor workers(nthreads);
   //Process the jobs and get a vector of the output
   std::vector<TH1F> process_output = workers.Map(workItem, input_arguments_multithreading);

   //THIS IS SO JANKY
   //I should consider posting on the CERN Root forums to get a better solution
   //in steps
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

   return h2;
}

DataStruct<TH1F>* score_lineal_histogram(TString filepath, float_t scoring_square_half_length, float_t scoring_sphere_spacing, float_t scoring_sphere_diameter, int track_oversample_factor, float_t CPE_range, Long_t random_seed = time(NULL),Int_t nhistoriestoanalyze = 0)
{
	/*
	DESCRIPTION OF THE ALGORITHM:
	1.) Rotation the edep point by a random rotation matrix.
	2.) Shift the edep point by a range
	3.) Oversample tracks by track_oversample_factor number of times
	4.) Do not try to place the track in a volume if the track is outside of the box region
	5.) Calculate lineal energy from that
	*/
	//ROOT::EnableImplicitMT(2); //amazingly, this slows things down
	//Initialize the local variables where the branches will point
	Double_t particle, edep;
	SVector3 particle_position = SVector3(); 
	int historyID;

	//Open the file, retrieve the tree, and determine the number of histories
	TFile *f = TFile::Open(filepath);
	TTree *microdosimetry;
	f->GetObject("microdosimetry",microdosimetry);
	microdosimetry->SetCacheSize(10000000U); //10? MB
	microdosimetry->SetBranchStatus("*",kFALSE);
	microdosimetry->AddBranchToCache("*",false);
	microdosimetry->AddBranchToCache("x",true);
	microdosimetry->AddBranchToCache("y",true);
	microdosimetry->AddBranchToCache("z",true);
	microdosimetry->AddBranchToCache("totalEnergyDeposit",true);
	microdosimetry->AddBranchToCache("eventID",true);
	microdosimetry->SetBranchStatus("x",kTRUE);
	microdosimetry->SetBranchStatus("y",kTRUE);
	microdosimetry->SetBranchStatus("z",kTRUE);
	microdosimetry->SetBranchStatus("totalEnergyDeposit",kTRUE);
	microdosimetry->SetBranchStatus("eventID",kTRUE);
	int total_histories = determine_number_of_events(microdosimetry); //have to do this here before the event ID branch pointer set below

	//Set up pointers from your tree branches to the local variables that you want to store them in
	microdosimetry->SetBranchAddress("x",&particle_position[0]);
	microdosimetry->SetBranchAddress("y",&particle_position[1]);
	microdosimetry->SetBranchAddress("z",&particle_position[2]);
	microdosimetry->SetBranchAddress("totalEnergyDeposit",&edep);
	microdosimetry->SetBranchAddress("eventID",&historyID);

	//Create the output struct
	DataStruct<TH1F>* output = new DataStruct<TH1F>(TH1F("Lineal energy histogram", "y*f(y)", 200, -2,1.5));

	//Initialize the geometry
	long long int index;
	int num_spheres_linear = TMath::Ceil(((scoring_square_half_length*2)/scoring_sphere_spacing)); //this is how many spheres there will be in a line
	if  (num_spheres_linear % 2 == 0) {cout << "even, adding 1 extra sphere" << endl; num_spheres_linear += 1;} //This is so there is always a central, 0,0,0 sphere
	long long int num_spheres_total = TMath::Power((num_spheres_linear),3);
	float_t top_sphere_offset = -(((float(num_spheres_linear))/2)-0.5)*scoring_sphere_spacing;//So the furthest sphere away in X,Y, or Z will be number of spheres plus the half center sphere away from the center
	std::unordered_map<int,double> energy_map; //define an unordered dictionary to hold edep values and associated volume

	//Initialize thie lineal energy histogram
	//output->data = TH1F("Lineal energy histogram", "y*f(y)", 200, -2,1.5); 
	BinLogX(&output->data); //transform the bins to logarithmic

	//Set values to initialize loop
	int current_history = 0;
	int entry_offset = 0; //the entry_offset
	int k = 0;
	int inside = 0; //a counter for the number of edeps within the box AND within a sphere
	int outside = 0; //a counter for the number of edeps within the box BUT outside of a sphere
	srand(random_seed); //Give current time as a seed to random number generator
	if (nhistoriestoanalyze == 0) {nhistoriestoanalyze = total_histories;} //if the default number of histories to analyze is unchanged, then analyze all
	microdosimetry->GetEntry(0); //Grab the first value in the tree


	#if VERBOSE == 1
	cout << "Number of spheres in a line: " << num_spheres_linear << endl;
	cout << "Total number of spheres in a cube: " << num_spheres_total << endl;
	cout << "largest sphere offset from central sphere: " << top_sphere_offset << endl;
	cout << "Total number of entries in file: " << microdosimetry->GetEntries() << endl;
	#endif

	//DESCRIPTION OF THE LOOP BELOW:
	//i increments once after a history (i.e. all the edeps associated with a track) has been itereated through oversample(j) # of times
	//j loops through the history oversample number of times
	//This is the most confusing part of the iterator,
	//k+entry_offset gives you the current edep point you're looking at
	//after a history has been iterated over, the entry_offset moves up to last value of k
	//where the last value of k is the first edep value of the NEXT history

	//POSSIBLE OPTIMIZATION:
	//I think that moving the j iteration to further out, so I only access the objects in memory once
	//and shift them repeatedly would be faster
	//It's actually not so simple: because I would have to hold all of the different shifts and rotation matrices in memory at the SAME time to do this
	//under the current method, I only calculate the shift once and rotation matrix once and use it repeatedly
	//Maybe I could make it such that I calculate a list of shifts and list of rotation matrices though, that could work and MAY be more efficient

	for(int i = 0; i<nhistoriestoanalyze;i++)
	{
		for (int j = 0; j<track_oversample_factor;j++) //this is an iterator over 1 event, oversample number of times
		{
			//Create new x,y,z shifts
			//These are constrained to be between -CPE_Range to + CPE_range
			SVector3 position_shifts = SVector3(((rand())*CPE_range*2/(RAND_MAX))-CPE_range,((rand())*CPE_range*2/(RAND_MAX))-CPE_range,((rand())*CPE_range*2/(RAND_MAX))-CPE_range);

			//Generate a new random rotation matrix
			SMatrix33 rotation_matrix = uniform_random_rotation_matrix(rand()/RAND_MAX,rand()/RAND_MAX,rand()/RAND_MAX);

			bool track_complete = false;
			k = 0;
			microdosimetry->GetEntry(entry_offset+k);

			while(track_complete == false) //this will keep iterating until the current track is complete
			{
				if (edep > 0) //this iterator actually is necessary, turns out there are some 0 edep events in the track file
				{
					//Transform particle posityion by the random rotation. Have to do this before shifting so the origin remains the same
					SVector3 position_rotated = rotation_matrix*particle_position;

					//Shift the rotated particle position by a random x,y,z shift
					SVector3 position_shifted = position_rotated+position_shifts;

					//Check if inside box
					if (abs(position_shifted[0]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2) && abs(position_shifted[1]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2) && abs(position_shifted[2]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2)) // if inside box
					{

						//Convert from x,y,z to index position
						SIntegerVector3 position_indices = SIntegerVector3(TMath::Nint((position_shifted[0]-top_sphere_offset)/scoring_sphere_spacing),TMath::Nint((position_shifted[1]-top_sphere_offset)/scoring_sphere_spacing),TMath::Nint((position_shifted[2]-top_sphere_offset)/scoring_sphere_spacing));

						//Figure out if x,y,z coordinate is within sphere
						//top_sphere_offset+(x_index*scoring_sphere_spacing) should give you the center of the sphere closest to your coordinate
						SVector3 position_difference_from_nearest_sphere = SVector3(position_indices[0],position_indices[1],position_indices[2]);
						position_difference_from_nearest_sphere = position_difference_from_nearest_sphere*scoring_sphere_spacing;
						position_difference_from_nearest_sphere = position_difference_from_nearest_sphere+top_sphere_offset-position_shifted;

						if(ROOT::Math::Mag(position_difference_from_nearest_sphere) <= (scoring_sphere_diameter/2))
						{
							//Okay you are inside the sphere
							index = position_indices[0] + (position_indices[1]*(num_spheres_linear)) + position_indices[2]*TMath::Power((num_spheres_linear),2); //Keep in mind that for the index it starts counting at zero
							double_t lineal_energy = edep/((2./3.)*scoring_sphere_diameter); //this should be ev/nm which is same a kev/um
							energy_map[index] += lineal_energy; //energy_map[index] += edep;
							inside += 1;

						} else {outside += 1;} 
					}
				}

				k++;
				microdosimetry->GetEntry(entry_offset+k); //increment to next interaction point

				if(historyID != current_history) //We have reached the end of the current oversampled track
				{
					track_complete = true;
				}

			}
			//Has finished iterating over current oversampled event. Output data
			std::unordered_map<int,double_t>::iterator it = energy_map.begin();
 			while (it != energy_map.cend())
 			{
 				output->data.Fill(it->second);
 				it++;
 			}
 			energy_map.clear(); //empty the map for next event
 		}
 		//Now we have completely finished looping over the oversampled event
 		current_history = historyID;
 		entry_offset += k;
	}
	//We have finished the analysis, close the file
	f->Close();

	#if VERBOSE == 1 
	cout << "EDeps in spheres: " << inside << " EDeps in box but outside spheres: " << outside << endl; //If you score outside/over inside should give: pi = 6/((outs/ints)+1)
	#endif

	PMF_to_PDF(&output->data);
	Prepare_for_Semilog(&output->data);

	//Fill the information strings
	output->information_strings["simulation_track_structure_filename"] = filepath;
	output->information_strings["simulation_random_seed"] = to_string(random_seed);
	output->information_strings["simulation_nhistories_analyzed"] = to_string(nhistoriestoanalyze);
	output->information_strings["simulation_nevents_analyzed"] = to_string(entry_offset);
	output->information_strings["simulation_noversamples"] = to_string(track_oversample_factor);
	output->information_strings["simulation_shift_length"] = to_string(CPE_range);
	output->information_strings["simulation_rotation_enabled"] = "True";
	output->information_strings["scoring_box_half_length"] = to_string(scoring_square_half_length);
	output->information_strings["scoring_sphere_diameter"] = to_string(scoring_sphere_diameter);
	output->information_strings["scoring_sphere_spacing"] = to_string(scoring_sphere_spacing);
	return output;
}
DataStruct<TH1F>* score_lineal_histogram_no_oversample(TString filepath, float_t scoring_square_half_length, float_t scoring_sphere_spacing, float_t scoring_sphere_diameter, int track_oversample_factor, float_t CPE_range, Long_t random_seed = time(NULL),Int_t nhistoriestoanalyze = 0)
{
	/*
	DESCRIPTION OF THE ALGORITHM:
	1.) Rotation the edep point by a random rotation matrix.
	2.) Shift the edep point by a range
	3.) Oversample tracks by track_oversample_factor number of times
	4.) Do not try to place the track in a volume if the track is outside of the box region
	5.) Calculate lineal energy from that
	*/
	//ROOT::EnableImplicitMT(2); //amazingly, this slows things down
	//Initialize the local variables where the branches will point
	Double_t particle, edep;
	SVector3 particle_position = SVector3(); 
	int historyID;

	//Open the file, retrieve the tree, and determine the number of histories
	TFile *f = TFile::Open(filepath);
	TTree *microdosimetry;
	f->GetObject("microdosimetry",microdosimetry);
	microdosimetry->SetCacheSize(10000000U); //10? MB
	microdosimetry->SetBranchStatus("*",kFALSE);
	microdosimetry->AddBranchToCache("*",false);
	microdosimetry->AddBranchToCache("x",true);
	microdosimetry->AddBranchToCache("y",true);
	microdosimetry->AddBranchToCache("z",true);
	microdosimetry->AddBranchToCache("totalEnergyDeposit",true);
	microdosimetry->AddBranchToCache("eventID",true);
	microdosimetry->SetBranchStatus("x",kTRUE);
	microdosimetry->SetBranchStatus("y",kTRUE);
	microdosimetry->SetBranchStatus("z",kTRUE);
	microdosimetry->SetBranchStatus("totalEnergyDeposit",kTRUE);
	microdosimetry->SetBranchStatus("eventID",kTRUE);
	int total_histories = determine_number_of_events(microdosimetry); //have to do this here before the event ID branch pointer set below

	//Set up pointers from your tree branches to the local variables that you want to store them in
	microdosimetry->SetBranchAddress("x",&particle_position[0]);
	microdosimetry->SetBranchAddress("y",&particle_position[1]);
	microdosimetry->SetBranchAddress("z",&particle_position[2]);
	microdosimetry->SetBranchAddress("totalEnergyDeposit",&edep);
	microdosimetry->SetBranchAddress("eventID",&historyID);

	//Create the output struct
	DataStruct<TH1F>* output = new DataStruct<TH1F>(TH1F("Lineal energy histogram", "y*f(y)", 200, -2,1.5));

	//Initialize the geometry
	long long int index;
	int num_spheres_linear = TMath::Ceil(((scoring_square_half_length*2)/scoring_sphere_spacing)); //this is how many spheres there will be in a line
	if  (num_spheres_linear % 2 == 0) {cout << "even, adding 1 extra sphere" << endl; num_spheres_linear += 1;} //This is so there is always a central, 0,0,0 sphere
	long long int num_spheres_total = TMath::Power((num_spheres_linear),3);
	float_t top_sphere_offset = -(((float(num_spheres_linear))/2)-0.5)*scoring_sphere_spacing;//So the furthest sphere away in X,Y, or Z will be number of spheres plus the half center sphere away from the center
	std::unordered_map<int,double> energy_map; //define an unordered dictionary to hold edep values and associated volume

	//Initialize thie lineal energy histogram
	//output->data = TH1F("Lineal energy histogram", "y*f(y)", 200, -2,1.5); 
	BinLogX(&output->data); //transform the bins to logarithmic

	//Set values to initialize loop
	int current_history = 0;
	int entry_offset = 0; //the entry_offset
	int k = 0;
	int inside = 0; //a counter for the number of edeps within the box AND within a sphere
	int outside = 0; //a counter for the number of edeps within the box BUT outside of a sphere
	srand(random_seed); //Give current time as a seed to random number generator
	if (nhistoriestoanalyze == 0) {nhistoriestoanalyze = total_histories;} //if the default number of histories to analyze is unchanged, then analyze all
	microdosimetry->GetEntry(0); //Grab the first value in the tree


	#if VERBOSE == 1
	cout << "Number of spheres in a line: " << num_spheres_linear << endl;
	cout << "Total number of spheres in a cube: " << num_spheres_total << endl;
	cout << "largest sphere offset from central sphere: " << top_sphere_offset << endl;
	cout << "Total number of entries in file: " << microdosimetry->GetEntries() << endl;
	#endif

	//DESCRIPTION OF THE LOOP BELOW:
	//i increments once after a history (i.e. all the edeps associated with a track) has been itereated through oversample(j) # of times
	//j loops through the history oversample number of times
	//This is the most confusing part of the iterator,
	//k+entry_offset gives you the current edep point you're looking at
	//after a history has been iterated over, the entry_offset moves up to last value of k
	//where the last value of k is the first edep value of the NEXT history

	//POSSIBLE OPTIMIZATION:
	//I think that moving the j iteration to further out, so I only access the objects in memory once
	//and shift them repeatedly would be faster
	//It's actually not so simple: because I would have to hold all of the different shifts and rotation matrices in memory at the SAME time to do this
	//under the current method, I only calculate the shift once and rotation matrix once and use it repeatedly
	//Maybe I could make it such that I calculate a list of shifts and list of rotation matrices though, that could work and MAY be more efficient

	for(int i = 0; i<nhistoriestoanalyze;i++)
	{
		
		//Create new x,y,z shifts
		//These are constrained to be between -CPE_Range to + CPE_range
		SVector3 position_shifts = SVector3(((rand())*CPE_range*2/(RAND_MAX))-CPE_range,((rand())*CPE_range*2/(RAND_MAX))-CPE_range,((rand())*CPE_range*2/(RAND_MAX))-CPE_range);

		//Generate a new random rotation matrix
		SMatrix33 rotation_matrix = uniform_random_rotation_matrix(rand()/RAND_MAX,rand()/RAND_MAX,rand()/RAND_MAX);

		bool track_complete = false;
		k = 0;
		microdosimetry->GetEntry(entry_offset+k);

		while(track_complete == false) //this will keep iterating until the current track is complete
		{
			if (edep > 0) //this iterator actually is necessary, turns out there are some 0 edep events in the track file
			{
				//Transform particle posityion by the random rotation. Have to do this before shifting so the origin remains the same
				SVector3 position_rotated = rotation_matrix*particle_position;

				//Shift the rotated particle position by a random x,y,z shift
				SVector3 position_shifted = position_rotated+position_shifts;

				//Check if inside box
				if (abs(position_shifted[0]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2) && abs(position_shifted[1]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2) && abs(position_shifted[2]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2)) // if inside box
				{

					//Convert from x,y,z to index position
					SIntegerVector3 position_indices = SIntegerVector3(TMath::Nint((position_shifted[0]-top_sphere_offset)/scoring_sphere_spacing),TMath::Nint((position_shifted[1]-top_sphere_offset)/scoring_sphere_spacing),TMath::Nint((position_shifted[2]-top_sphere_offset)/scoring_sphere_spacing));

					//Figure out if x,y,z coordinate is within sphere
					//top_sphere_offset+(x_index*scoring_sphere_spacing) should give you the center of the sphere closest to your coordinate
					SVector3 position_difference_from_nearest_sphere = SVector3(position_indices[0],position_indices[1],position_indices[2]);
					position_difference_from_nearest_sphere = position_difference_from_nearest_sphere*scoring_sphere_spacing;
					position_difference_from_nearest_sphere = position_difference_from_nearest_sphere+top_sphere_offset-position_shifted;

					if(ROOT::Math::Mag(position_difference_from_nearest_sphere) <= (scoring_sphere_diameter/2))
					{
						//Okay you are inside the sphere
						index = position_indices[0] + (position_indices[1]*(num_spheres_linear)) + position_indices[2]*TMath::Power((num_spheres_linear),2); //Keep in mind that for the index it starts counting at zero
						double_t lineal_energy = edep/((2./3.)*scoring_sphere_diameter); //this should be ev/nm which is same a kev/um
						energy_map[index] += lineal_energy; //energy_map[index] += edep;
						inside += 1;

					} else {outside += 1;} 
				}
			}

			k++;
			microdosimetry->GetEntry(entry_offset+k); //increment to next interaction point

			if(historyID != current_history) //We have reached the end of the current oversampled track
			{
				track_complete = true;
			}

		}
		//Has finished iterating over current oversampled event. Output data
		std::unordered_map<int,double_t>::iterator it = energy_map.begin();
			while (it != energy_map.cend())
			{
				output->data.Fill(it->second);
				it++;
			}
			energy_map.clear(); //empty the map for next event
 		
 		//Now we have completely finished looping over the oversampled event
 		current_history = historyID;
 		entry_offset += k;
	}
	//We have finished the analysis, close the file
	f->Close();

	#if VERBOSE == 1 
	cout << "EDeps in spheres: " << inside << " EDeps in box but outside spheres: " << outside << endl; //If you score outside/over inside should give: pi = 6/((outs/ints)+1)
	#endif

	PMF_to_PDF(&output->data);
	Prepare_for_Semilog(&output->data);

	//Fill the information strings
	output->information_strings["simulation_track_structure_filename"] = filepath;
	output->information_strings["simulation_random_seed"] = to_string(random_seed);
	output->information_strings["simulation_nhistories_analyzed"] = to_string(nhistoriestoanalyze);
	output->information_strings["simulation_nevents_analyzed"] = to_string(entry_offset);
	output->information_strings["simulation_noversamples"] = to_string(track_oversample_factor);
	output->information_strings["simulation_shift_length"] = to_string(CPE_range);
	output->information_strings["simulation_rotation_enabled"] = "True";
	output->information_strings["scoring_box_half_length"] = to_string(scoring_square_half_length);
	output->information_strings["scoring_sphere_diameter"] = to_string(scoring_sphere_diameter);
	output->information_strings["scoring_sphere_spacing"] = to_string(scoring_sphere_spacing);
	return output;
}

DataStruct<std::vector<double>>* score_lineal_list(TString filepath, float_t scoring_square_half_length, float_t scoring_sphere_spacing, float_t scoring_sphere_diameter, int track_oversample_factor, float_t CPE_range, Long_t random_seed = time(NULL),Int_t nhistoriestoanalyze = 0)
{
	/*
	DESCRIPTION OF THE ALGORITHM:
	1.) Rotation the edep point by a random rotation matrix.
	2.) Shift the edep point by a range
	3.) Oversample tracks by track_oversample_factor number of times
	4.) Do not try to place the track in a volume if the track is outside of the box region
	5.) Calculate lineal energy from that
	*/

	//Initialize the local variables where the branches will point
	Double_t particle, edep;
	SVector3 particle_position = SVector3(); 
	int historyID;

	//Open the file, retrieve the tree, and determine the number of histories
	TFile *f = TFile::Open(filepath);
	TTree *microdosimetry;
	f->GetObject("microdosimetry",microdosimetry);
	int total_histories = determine_number_of_events(microdosimetry); //have to do this here before the event ID branch pointer set below

	//Set up pointers from your tree branches to the local variables that you want to store them in
	microdosimetry->SetBranchAddress("x",&particle_position[0]);
	microdosimetry->SetBranchAddress("y",&particle_position[1]);
	microdosimetry->SetBranchAddress("z",&particle_position[2]);
	microdosimetry->SetBranchAddress("totalEnergyDeposit",&edep);
	microdosimetry->SetBranchAddress("eventID",&historyID);

	//Create the output struct
	DataStruct<std::vector<double>>* output = new DataStruct<std::vector<double>>(std::vector<double>());

	//Initialize the geometry
	long long int index;
	int num_spheres_linear = TMath::Ceil(((scoring_square_half_length*2)/scoring_sphere_spacing)); //this is how many spheres there will be in a line
	if  (num_spheres_linear % 2 == 0) {cout << "even, adding 1 extra sphere" << endl; num_spheres_linear += 1;} //This is so there is always a central, 0,0,0 sphere
	long long int num_spheres_total = TMath::Power((num_spheres_linear),3);
	float_t top_sphere_offset = -(((float(num_spheres_linear))/2)-0.5)*scoring_sphere_spacing;//So the furthest sphere away in X,Y, or Z will be number of spheres plus the half center sphere away from the center
	std::unordered_map<int,double> energy_map; //define an unordered dictionary to hold edep values and associated volume

	//Set values to initialize loop
	int current_history = 0;
	int entry_offset = 0; //the entry_offset
	int k = 0;
	int inside = 0; //a counter for the number of edeps within the box AND within a sphere
	int outside = 0; //a counter for the number of edeps within the box BUT outside of a sphere
	srand(random_seed); //Give current time as a seed to random number generator
	if (nhistoriestoanalyze == 0) {nhistoriestoanalyze = total_histories;} //if the default number of histories to analyze is unchanged, then analyze all
	microdosimetry->GetEntry(0); //Grab the first value in the tree


	#if VERBOSE == 1
	cout << "Number of spheres in a line: " << num_spheres_linear << endl;
	cout << "Total number of spheres in a cube: " << num_spheres_total << endl;
	cout << "largest sphere offset from central sphere: " << top_sphere_offset << endl;
	cout << "Total number of entries in file: " << microdosimetry->GetEntries() << endl;
	#endif

	//DESCRIPTION OF THE LOOP BELOW:
	//i increments once after a history (i.e. all the edeps associated with a track) has been itereated through oversample(j) # of times
	//j loops through the history oversample number of times
	//This is the most confusing part of the iterator,
	//k+entry_offset gives you the current edep point you're looking at
	//after a history has been iterated over, the entry_offset moves up to last value of k
	//where the last value of k is the first edep value of the NEXT history

	//POSSIBLE OPTIMIZATION:
	//I think that moving the j iteration to further out, so I only access the objects in memory once
	//and shift them repeatedly would be faster
	//It's actually not so simple: because I would have to hold all of the different shifts and rotation matrices in memory at the SAME time to do this
	//under the current method, I only calculate the shift once and rotation matrix once and use it repeatedly
	//Maybe I could make it such that I calculate a list of shifts and list of rotation matrices though, that could work and MAY be more efficient

	for(int i = 0; i<nhistoriestoanalyze;i++)
	{
		for (int j = 0; j<track_oversample_factor;j++) //this is an iterator over 1 event, oversample number of times
		{
			//Create new x,y,z shifts
			//These are constrained to be between -CPE_Range to + CPE_range
			SVector3 position_shifts = SVector3(((rand())*CPE_range*2/(RAND_MAX))-CPE_range,((rand())*CPE_range*2/(RAND_MAX))-CPE_range,((rand())*CPE_range*2/(RAND_MAX))-CPE_range);

			//Generate a new random rotation matrix
			SMatrix33 rotation_matrix = uniform_random_rotation_matrix(rand()/RAND_MAX,rand()/RAND_MAX,rand()/RAND_MAX);

			bool track_complete = false;
			k = 0;
			microdosimetry->GetEntry(entry_offset+k);

			while(track_complete == false) //this will keep iterating until the current track is complete
			{
				if (edep > 0) //this iterator actually is necessary, turns out there are some 0 edep events in the track file
				{
					//Transform particle posityion by the random rotation. Have to do this before shifting so the origin remains the same
					SVector3 position_rotated = rotation_matrix*particle_position;

					//Shift the rotated particle position by a random x,y,z shift
					SVector3 position_shifted = position_rotated+position_shifts;

					//Check if inside box
					if (abs(position_shifted[0]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2) && abs(position_shifted[1]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2) && abs(position_shifted[2]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2)) // if inside box
					{

						//Convert from x,y,z to index position
						SIntegerVector3 position_indices = SIntegerVector3(TMath::Nint((position_shifted[0]-top_sphere_offset)/scoring_sphere_spacing),TMath::Nint((position_shifted[1]-top_sphere_offset)/scoring_sphere_spacing),TMath::Nint((position_shifted[2]-top_sphere_offset)/scoring_sphere_spacing));

						//Figure out if x,y,z coordinate is within sphere
						//top_sphere_offset+(x_index*scoring_sphere_spacing) should give you the center of the sphere closest to your coordinate
						SVector3 position_difference_from_nearest_sphere = SVector3(position_indices[0],position_indices[1],position_indices[2]);
						position_difference_from_nearest_sphere = position_difference_from_nearest_sphere*scoring_sphere_spacing;
						position_difference_from_nearest_sphere = position_difference_from_nearest_sphere+top_sphere_offset-position_shifted;

						if(ROOT::Math::Mag(position_difference_from_nearest_sphere) <= (scoring_sphere_diameter/2))
						{
							//Okay you are inside the sphere
							index = position_indices[0] + (position_indices[1]*(num_spheres_linear)) + position_indices[2]*TMath::Power((num_spheres_linear),2); //Keep in mind that for the index it starts counting at zero
							double_t lineal_energy = edep/((2./3.)*scoring_sphere_diameter); //this should be ev/nm which is same a kev/um
							energy_map[index] += lineal_energy; //energy_map[index] += edep;
							inside += 1;

						} else {outside += 1;} 
					}
				}

				k++;
				microdosimetry->GetEntry(entry_offset+k); //increment to next interaction point

				if(historyID != current_history) //We have reached the end of the current oversampled track
				{
					track_complete = true;
				}

			}
			//Has finished iterating over current oversampled event. Output data
			std::unordered_map<int,double_t>::iterator it = energy_map.begin();
 			while (it != energy_map.cend())
 			{
 				output->data.push_back(it->second);
 				it++;
 			}
 			energy_map.clear(); //empty the map for next event
 		}
 		//Now we have completely finished looping over the oversampled event
 		current_history = historyID;
 		entry_offset += k;
	}
	//We have finished the analysis, close the file
	f->Close();

	#if VERBOSE == 1 
	cout << "EDeps in spheres: " << inside << " EDeps in box but outside spheres: " << outside << endl; //If you score outside/over inside should give: pi = 6/((outs/ints)+1)
	#endif

	//Fill the information strings
	output->information_strings["simulation_track_structure_filename"] = filepath;
	output->information_strings["simulation_random_seed"] = to_string(random_seed);
	output->information_strings["simulation_nhistories_analyzed"] = to_string(nhistoriestoanalyze);
	output->information_strings["simulation_nevents_analyzed"] = to_string(entry_offset);
	output->information_strings["simulation_noversamples"] = to_string(track_oversample_factor);
	output->information_strings["simulation_shift_length"] = to_string(CPE_range);
	output->information_strings["simulation_rotation_enabled"] = "True";
	output->information_strings["scoring_box_half_length"] = to_string(scoring_square_half_length);
	output->information_strings["scoring_sphere_diameter"] = to_string(scoring_sphere_diameter);
	output->information_strings["scoring_sphere_spacing"] = to_string(scoring_sphere_spacing);
	return output;
}

/*DataStruct<std::unordered_map<int,std::vector<double>>>* score_lineal_paired_list(TString filepath, float_t scoring_square_half_length, float_t scoring_sphere_spacing, float_t scoring_sphere_diameter, int track_oversample_factor, float_t CPE_range, Long_t random_seed = time(NULL),Int_t nhistoriestoanalyze = 0)
{
	
	//Initialize the local variables where the branches will point
	Double_t particle, edep;
	SVector3 particle_position = SVector3(); 
	int historyID;

	//Open the file, retrieve the tree, and determine the number of histories
	TFile *f = TFile::Open(filepath);
	TTree *microdosimetry;
	f->GetObject("microdosimetry",microdosimetry);
	int total_histories = determine_number_of_events(microdosimetry); //have to do this here before the event ID branch pointer set below

	//Set up pointers from your tree branches to the local variables that you want to store them in
	microdosimetry->SetBranchAddress("x",&particle_position[0]);
	microdosimetry->SetBranchAddress("y",&particle_position[1]);
	microdosimetry->SetBranchAddress("z",&particle_position[2]);
	microdosimetry->SetBranchAddress("totalEnergyDeposit",&edep);
	microdosimetry->SetBranchAddress("eventID",&historyID);

	//Create the output struct
	DataStruct<std::unordered_map<int,std::vector<double>>>* output = new DataStruct<std::unordered_map<int,std::vector<double>>>(std::unordered_map<int,std::vector<double>>());

	//Initialize the geometry
	long long int index;
	int num_spheres_linear = TMath::Ceil(((scoring_square_half_length*2)/scoring_sphere_spacing)); //this is how many spheres there will be in a line
	if  (num_spheres_linear % 2 == 0) {cout << "even, adding 1 extra sphere" << endl; num_spheres_linear += 1;} //This is so there is always a central, 0,0,0 sphere
	long long int num_spheres_total = TMath::Power((num_spheres_linear),3);
	float_t top_sphere_offset = -(((float(num_spheres_linear))/2)-0.5)*scoring_sphere_spacing;//So the furthest sphere away in X,Y, or Z will be number of spheres plus the half center sphere away from the center
	std::unordered_map<int,double> energy_map; //define an unordered dictionary to hold edep values and associated volume

	//Set values to initialize loop
	int current_history = 0;
	int entry_offset = 0; //the entry_offset
	int k = 0;
	int inside = 0; //a counter for the number of edeps within the box AND within a sphere
	int outside = 0; //a counter for the number of edeps within the box BUT outside of a sphere
	srand(random_seed); //Give current time as a seed to random number generator
	if (nhistoriestoanalyze == 0) {nhistoriestoanalyze = total_histories;} //if the default number of histories to analyze is unchanged, then analyze all
	microdosimetry->GetEntry(0); //Grab the first value in the tree


	#if VERBOSE == 1
	cout << "Number of spheres in a line: " << num_spheres_linear << endl;
	cout << "Total number of spheres in a cube: " << num_spheres_total << endl;
	cout << "largest sphere offset from central sphere: " << top_sphere_offset << endl;
	cout << "Total number of entries in file: " << microdosimetry->GetEntries() << endl;
	#endif

	//DESCRIPTION OF THE LOOP BELOW:
	//i increments once after a history (i.e. all the edeps associated with a track) has been itereated through oversample(j) # of times
	//j loops through the history oversample number of times
	//This is the most confusing part of the iterator,
	//k+entry_offset gives you the current edep point you're looking at
	//after a history has been iterated over, the entry_offset moves up to last value of k
	//where the last value of k is the first edep value of the NEXT history

	//POSSIBLE OPTIMIZATION:
	//I think that moving the j iteration to further out, so I only access the objects in memory once
	//and shift them repeatedly would be faster
	//It's actually not so simple: because I would have to hold all of the different shifts and rotation matrices in memory at the SAME time to do this
	//under the current method, I only calculate the shift once and rotation matrix once and use it repeatedly
	//Maybe I could make it such that I calculate a list of shifts and list of rotation matrices though, that could work and MAY be more efficient

	for(int i = 0; i<nhistoriestoanalyze;i++)
	{
		for (int j = 0; j<track_oversample_factor;j++) //this is an iterator over 1 event, oversample number of times
		{
			//Create new x,y,z shifts
			//These are constrained to be between -CPE_Range to + CPE_range
			SVector3 position_shifts = SVector3(((rand())*CPE_range*2/(RAND_MAX))-CPE_range,((rand())*CPE_range*2/(RAND_MAX))-CPE_range,((rand())*CPE_range*2/(RAND_MAX))-CPE_range);

			//Generate a new random rotation matrix
			SMatrix33 rotation_matrix = uniform_random_rotation_matrix(rand()/RAND_MAX,rand()/RAND_MAX,rand()/RAND_MAX);

			bool track_complete = false;
			k = 0;
			microdosimetry->GetEntry(entry_offset+k);

			while(track_complete == false) //this will keep iterating until the current track is complete
			{
				if (edep > 0) //this iterator actually is necessary, turns out there are some 0 edep events in the track file
				{
					//Transform particle posityion by the random rotation. Have to do this before shifting so the origin remains the same
					SVector3 position_rotated = rotation_matrix*particle_position;

					//Shift the rotated particle position by a random x,y,z shift
					SVector3 position_shifted = position_rotated+position_shifts;

					//Check if inside box
					if (abs(position_shifted[0]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2) && abs(position_shifted[1]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2) && abs(position_shifted[2]) < abs(top_sphere_offset)+(scoring_sphere_diameter/2)) // if inside box
					{

						//Convert from x,y,z to index position
						SIntegerVector3 position_indices = SIntegerVector3(TMath::Nint((position_shifted[0]-top_sphere_offset)/scoring_sphere_spacing),TMath::Nint((position_shifted[1]-top_sphere_offset)/scoring_sphere_spacing),TMath::Nint((position_shifted[2]-top_sphere_offset)/scoring_sphere_spacing));

						//Figure out if x,y,z coordinate is within sphere
						//top_sphere_offset+(x_index*scoring_sphere_spacing) should give you the center of the sphere closest to your coordinate
						SVector3 position_difference_from_nearest_sphere = SVector3(position_indices[0],position_indices[1],position_indices[2]);
						position_difference_from_nearest_sphere = position_difference_from_nearest_sphere*scoring_sphere_spacing;
						position_difference_from_nearest_sphere = position_difference_from_nearest_sphere+top_sphere_offset-position_shifted;

						if(ROOT::Math::Mag(position_difference_from_nearest_sphere) <= (scoring_sphere_diameter/2))
						{
							//Okay you are inside the sphere
							index = position_indices[0] + (position_indices[1]*(num_spheres_linear)) + position_indices[2]*TMath::Power((num_spheres_linear),2); //Keep in mind that for the index it starts counting at zero
							double_t lineal_energy = edep/((2./3.)*scoring_sphere_diameter); //this should be ev/nm which is same a kev/um
							energy_map[index] += lineal_energy; //energy_map[index] += edep;
							inside += 1;

						} else {outside += 1;} 
					}
				}

				k++;
				microdosimetry->GetEntry(entry_offset+k); //increment to next interaction point

				if(historyID != current_history) //We have reached the end of the current oversampled track
				{
					track_complete = true;
				}

			}
			//Has finished iterating over current oversampled event. Output data
			std::unordered_map<int,double_t>::iterator it = energy_map.begin();
 			while (it != energy_map.cend())
 			{
 				output->data[index].push_back(it->second);
 				it++;
 			}
 			energy_map.clear(); //empty the map for next event
 		}
 		//Now we have completely finished looping over the oversampled event
 		current_history = historyID;
 		entry_offset += k;
	}
	//We have finished the analysis, close the file
	f->Close();

	#if VERBOSE == 1 
	cout << "EDeps in spheres: " << inside << " EDeps in box but outside spheres: " << outside << endl; //If you score outside/over inside should give: pi = 6/((outs/ints)+1)
	#endif

	//Fill the information strings
	output->information_strings["simulation_track_structure_filename"] = filepath;
	output->information_strings["simulation_random_seed"] = to_string(random_seed);
	output->information_strings["simulation_nhistories_analyzed"] = to_string(nhistoriestoanalyze);
	output->information_strings["simulation_nevents_analyzed"] = to_string(entry_offset);
	output->information_strings["simulation_noversamples"] = to_string(track_oversample_factor);
	output->information_strings["simulation_shift_length"] = to_string(CPE_range);
	output->information_strings["simulation_rotation_enabled"] = "True";
	output->information_strings["scoring_box_half_length"] = to_string(scoring_square_half_length);
	output->information_strings["scoring_sphere_diameter"] = to_string(scoring_sphere_diameter);
	output->information_strings["scoring_sphere_spacing"] = to_string(scoring_sphere_spacing);
	return output;
}*/

void analyze_lineal()
{
	int start_time = time(0); cout << "ROOT Program Beginning" << endl; 	

	//Recall that sizes are in nm, so e3=um,e6=mm,e7=cm
	// 1.) Scoring square half length 2.) Sphere spacing 3.) Sphere diameter 4.) Track oversample factor 5.) CPE range
	DataStruct<TH1F>* output = score_lineal_histogram_no_oversample("/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm.root",5e6, 5e3, 5e3, 5, 50e3, time(NULL), 10000); //7.4 um is the mean nucleus diameter for the malignant nuclei of ptn. 1 slide 1 (i.e. it should match the results shown in the manuscript)
	output->save("/home/joseph/Documents/M1_local/mar21");

	//Plotting
	THStack *histo_stack = new THStack("histograms","");
	histo_stack->Add(&output->data);
	histo_stack->Draw("nostack"); //Draw histogram
	gPad->SetLogx(); //Set the logarithmic axes appropriately

	//Remember if you're making a bigger program you have to do memory management
	//delete histo_stack; //delete output; /Don't do these if you want to plot the histogram otherwise it won't appear
	
	int end_time = time(0); cout << "ROOT Program Ending. Seconds elapsed: " << (end_time-start_time) << endl;	
}

void analyze_lineal_list()
{
	int start_time = time(0); cout << "ROOT Program Beginning" << endl; 	

	//Recall that sizes are in nm, so e3=um,e6=mm,e7=cm
	// 1.) Scoring square half length 2.) Sphere spacing 3.) Sphere diameter 4.) Track oversample factor 5.) CPE range
	DataStruct<std::vector<double>>* output = score_lineal_list("/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm.root",5e6, 5e3, 5e3, 5, 50e3, time(NULL), 100); //7.4 um is the mean nucleus diameter for the malignant nuclei of ptn. 1 slide 1 (i.e. it should match the results shown in the manuscript)
	output->save("/home/joseph/Documents/M1_local/mar21");

	//Plotting

	//Remember if you're making a bigger program you have to do memory management
	//delete output; /Don't do these if you want to plot the histogram otherwise it won't appear
	
	int end_time = time(0); cout << "ROOT Program Ending. Seconds elapsed: " << (end_time-start_time) << endl;	
}

/*void analyze_lineal_multithreaded()
{
	int start_time = time(0); cout << "ROOT Program Beginning" << endl; 	

	//Recall that sizes are in nm, so e3=um,e6=mm,e7=cm
	// 1.) Scoring square half length 2.) Sphere spacing 3.) Sphere diameter 4.) Track oversample factor 5.) CPE range
	score_lineal_histogram_multithreaded("/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm.root",5e6, 10e3, 5e3, 5, 50e3, time(NULL), 100); //7.4 um is the mean nucleus diameter for the malignant nuclei of ptn. 1 slide 1 (i.e. it should match the results shown in the manuscript)

	//Remember if you're making a bigger program you have to do memory management
	//delete histo_stack; //delete output; /Don't do these if you want to plot the histogram otherwise it won't appear
	
	int end_time = time(0); cout << "ROOT Program Ending. Seconds elapsed: " << (end_time-start_time) << endl;	
}*/

void analyze_lineal_multithreaded_explicit()
{
	int start_time = time(0); cout << "ROOT Program Beginning" << endl; 	

	//Recall that sizes are in nm, so e3=um,e6=mm,e7=cm
	// 1.) Scoring square half length 2.) Sphere spacing 3.) Sphere diameter 4.) CPE range
	//auto output = score_lineal_histogram_multithreaded_explicit("/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm_4vectortransformed_1000000.root",5e6, 5e3, 5e3, 50e3, 1, time(NULL), 1000000); //7.4 um is the mean nucleus diameter for the malignant nuclei of ptn. 1 slide 1 (i.e. it should match the results shown in the manuscript)
	//auto output = new TH1F(score_lineal_histogram_multithreaded_explicit("/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm_4vectortransformed_1000000.root",1e6, 8.82e3, 7.4e3, 5e6, 4, time(NULL), 10000));
	TObjArray output_array(0);
	Int_t noversamples = 1;
	//Plotting
	THStack *histo_stack = new THStack("histograms","");
	//IF YOU WANT TO HAVE GARBAGE COLLECTION. Set the list as the owner. Except that being the case, then your objects won't appear since they die and go out of scope.
	//output_array.SetOwner(true);
	for (int i = 0; i < noversamples; i++)
	{
		//output_array.Add(new TH1F(score_lineal_histogram_multithreaded_explicit("/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm_4vectortransformed_1000000.root",125e3, 17.5e3, 7.4e3, 5e6, 4, time(NULL), 1000000)));
		output_array.Add(new TH1F(score_lineal_histogram_multithreaded_explicit("/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/software/MicroTrackGenerator/output/co60_1cm_livermore_4vectortransformed_1000.root",5e6, 12e3, 12e3, 12e3, 4, time(NULL))));
	}
	TH1F* histo = (TH1F*)output_array[0]->Clone(); //I think the way clone works, is if the underlying object is on the heap, so is the clone. ALSO the original object takes ownership of the object. So if I delete the original object I lose the clone.
	for (int i = 1; i < noversamples; i++)
	{
		histo->Add((TH1F*)output_array[i]);
	}
	//TH1F* histo = (TH1F*)output_list[0].Clone();
	//histo->Merge(output_list);
	//TH1F* histo_heap = new TH1F;
	//histo_heap = histo;
	histo_stack->Draw(); //Draw histogram
	histo_stack->Add(histo);
	histo_stack->GetXaxis()->SetTitle("y [keV/um]"); 
	//histo_stack->GetXaxis()->SetTitle("edep_1 [keV]"); 
	histo_stack->GetXaxis()->CenterTitle();
	histo_stack->GetYaxis()->SetTitle("y*f(y)"); 
	//histo_stack->GetXaxis()->SetTitle("edep_1 * f(edep_1)"); 
	histo_stack->GetYaxis()->CenterTitle();
	gPad->Modified();
	gPad->SetLogx(); //Set the logarithmic axes appropriately

	//Remember if you're making a bigger program you have to do memory management
	//delete histo_stack; //delete output; /Don't do these if you want to plot the histogram otherwise it won't appear
	
	int end_time = time(0); cout << "ROOT Program Ending. Seconds elapsed: " << (end_time-start_time) << endl;	
}

void debug_verify_position_shifts()
{
	int start_time = time(0); cout << "ROOT Program Beginning" << endl; 	

	//Recall that sizes are in nm, so e3=um,e6=mm,e7=cm
	// 1.) Scoring square half length 2.) Sphere spacing 3.) Sphere diameter 4.) Track oversample factor 5.) CPE range
	//auto output = score_lineal_histogram_multithreaded_explicit("/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm_4vectortransformed_1000000.root",5e6, 5e3, 5e3, 50e3, 1, time(NULL), 1000000); //7.4 um is the mean nucleus diameter for the malignant nuclei of ptn. 1 slide 1 (i.e. it should match the results shown in the manuscript)
	//auto output = new TH1F(score_lineal_histogram_multithreaded_explicit("/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm_4vectortransformed_1000000.root",1e6, 8.82e3, 7.4e3, 5e6, 4, time(NULL), 10000));
	TObjArray output_array(0);
	Int_t noversamples = 1;
	//Plotting
	THStack *histo_stack = new THStack("histograms","");
	//IF YOU WANT TO HAVE GARBAGE COLLECTION. Set the list as the owner. Except that being the case, then your objects won't appear since they die and go out of scope.
	//output_array.SetOwner(true);
	for (int i = 0; i < noversamples; i++)
	{
		//output_array.Add(new TH1F(score_lineal_histogram_multithreaded_explicit("/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm_4vectortransformed_1000000.root",125e3, 17.5e3, 7.4e3, 5e6, 4, time(NULL), 1000000)));
		output_array.Add(new TH1F(debug_code_verify_position_shifts("/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm_4vectortransformed_1000000.root",5e6, 17.5e3, 7.4e3, 5e6, 4, time(NULL), 1000000)));
	}
	TH1F* histo = (TH1F*)output_array[0]->Clone(); //I think the way clone works, is if the underlying object is on the heap, so is the clone. ALSO the original object takes ownership of the object. So if I delete the original object I lose the clone.
	for (int i = 1; i < noversamples; i++)
	{
		histo->Add((TH1F*)output_array[i]);
	}
	//TH1F* histo = (TH1F*)output_list[0].Clone();
	//histo->Merge(output_list);
	//TH1F* histo_heap = new TH1F;
	//histo_heap = histo;
	histo_stack->Draw(); //Draw histogram
	histo_stack->Add(histo);
	histo_stack->GetXaxis()->SetTitle("Spatial Coordinate"); 
	histo_stack->GetXaxis()->CenterTitle();
	histo_stack->GetYaxis()->SetTitle("Frequency"); 
	histo_stack->GetYaxis()->CenterTitle();
	gPad->Modified();
	//gPad->SetLogx(); //Set the logarithmic axes appropriately

	//Remember if you're making a bigger program you have to do memory management
	//delete histo_stack; //delete output; /Don't do these if you want to plot the histogram otherwise it won't appear
	
	int end_time = time(0); cout << "ROOT Program Ending. Seconds elapsed: " << (end_time-start_time) << endl;	
}

void debug_verify_rotation_shifts()
{
	int start_time = time(0); cout << "ROOT Program Beginning" << endl; 	

	//Recall that sizes are in nm, so e3=um,e6=mm,e7=cm
	// 1.) Scoring square half length 2.) Sphere spacing 3.) Sphere diameter 4.) Track oversample factor 5.) CPE range
	//auto output = score_lineal_histogram_multithreaded_explicit("/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm_4vectortransformed_1000000.root",5e6, 5e3, 5e3, 50e3, 1, time(NULL), 1000000); //7.4 um is the mean nucleus diameter for the malignant nuclei of ptn. 1 slide 1 (i.e. it should match the results shown in the manuscript)
	//auto output = new TH1F(score_lineal_histogram_multithreaded_explicit("/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm_4vectortransformed_1000000.root",1e6, 8.82e3, 7.4e3, 5e6, 4, time(NULL), 10000));
	TObjArray output_array(0);
	Int_t noversamples = 1;
	//Plotting
	THStack *histo_stack = new THStack("histograms","");
	//IF YOU WANT TO HAVE GARBAGE COLLECTION. Set the list as the owner. Except that being the case, then your objects won't appear since they die and go out of scope.
	//output_array.SetOwner(true);
	for (int i = 0; i < noversamples; i++)
	{
		//output_array.Add(new TH1F(score_lineal_histogram_multithreaded_explicit("/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm_4vectortransformed_1000000.root",125e3, 17.5e3, 7.4e3, 5e6, 4, time(NULL), 1000000)));
		output_array.Add(new TH1F(debug_code_verify_rotation_shifts("/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm_4vectortransformed_1000000.root",5e6, 17.5e3, 7.4e3, 5e6, 4, time(NULL), 4,2)));
	}
	TH1F* histo = (TH1F*)output_array[0]->Clone(); //I think the way clone works, is if the underlying object is on the heap, so is the clone. ALSO the original object takes ownership of the object. So if I delete the original object I lose the clone.
	for (int i = 1; i < noversamples; i++)
	{
		histo->Add((TH1F*)output_array[i]);
	}
	//TH1F* histo = (TH1F*)output_list[0].Clone();
	//histo->Merge(output_list);
	//TH1F* histo_heap = new TH1F;
	//histo_heap = histo;
	histo_stack->Draw(); //Draw histogram
	histo_stack->Add(histo);
	histo_stack->GetXaxis()->SetTitle("Spatial Coordinate"); 
	histo_stack->GetXaxis()->CenterTitle();
	histo_stack->GetYaxis()->SetTitle("Frequency"); 
	histo_stack->GetYaxis()->CenterTitle();
	gPad->Modified();
	//gPad->SetLogx(); //Set the logarithmic axes appropriately

	//Remember if you're making a bigger program you have to do memory management
	//delete histo_stack; //delete output; /Don't do these if you want to plot the histogram otherwise it won't appear
	
	int end_time = time(0); cout << "ROOT Program Ending. Seconds elapsed: " << (end_time-start_time) << endl;	
}


void transform_file_to_4vector()
{
	TString filepath = "/home/joseph/Dropbox/Documents/Work/Projects/Alpha_Microdosimetry/software/MicroTrackGenerator/output/co60_1cm_livermore.root";
	transform_microdosimetry_file_to_4_vector(filepath,false,1000);
	/*
	TString filepath = "/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm_4vectortransformed_10000.root"
	TFile *f = TFile::Open(filepath);
	TTree *microdosimetry;
	f->GetObject("microdosimetry",microdosimetry);
	microdosimetry->Print();
	*/
}
__global__ void GPU_with_ROOT_test()
{
	printf("Yolo \n");
}

void superimpose_tracks()
{
	//analyze_lineal_multithreaded_explicit();
	cout << "Program continues to superimpose_tracks()" << endl;
	GPU_with_ROOT_test<<<1,12>>>();
	cudaDeviceSynchronize();
}
# ifndef __CINT__
int main()
{
	cout << "Program begins in main" << endl;
  superimpose_tracks();
  return 0;
}
# endif
/*void analyze_lineal_paired_list()
{
	int start_time = time(0); cout << "ROOT Program Beginning" << endl; 	

	//Recall that sizes are in nm, so e3=um,e6=mm,e7=cm
	// 1.) Scoring square half length 2.) Sphere spacing 3.) Sphere diameter 4.) Track oversample factor 5.) CPE range
	DataStruct<std::unordered_map<int,std::vector<double>>>* output = score_lineal_paired_list("/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm.root",5e6, 5e3, 5e3, 5, 50e3, time(NULL), 100); //7.4 um is the mean nucleus diameter for the malignant nuclei of ptn. 1 slide 1 (i.e. it should match the results shown in the manuscript)
	output->save("/home/joseph/Documents/M1_local/mar21");

	//Plotting

	//Remember if you're making a bigger program you have to do memory management
	//delete output; /Don't do these if you want to plot the histogram otherwise it won't appear
	
	int end_time = time(0); cout << "ROOT Program Ending. Seconds elapsed: " << (end_time-start_time) << endl;	
}*/

/*
Pseudo-code for multithreading:

Think this is the example I want to follow
https://root.cern.ch/doc/v612/mt103__fillNtupleFromMultipleThreads_8C.html

Interestingly, I wonder if conventional C++  threads will work for this purpose
or if I have to use built in root tools.

*/

/*
Scratch board
TFile *testfile = TFile::Open("/home/joseph/Documents/M1_local/mar21/co60_1cm_1617146107_100x5.root");
testfile->ls(); //to list everything inside

----
//time to add some svectors and smatrices
typedef ROOT::Math::SVector<Double_t,3> SVector3;
TString filepath = "/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm.root";	
TFile *f = TFile::Open(filepath);
TTree *microdosimetry;
f->GetObject("microdosimetry",microdosimetry);


SVector3 particle_location = SVector3();
microdosimetry->SetBranchAddress("x",&particle_location[0]);
microdosimetry->SetBranchAddress("y",&particle_location[1]);
microdosimetry->SetBranchAddress("z",&particle_location[2]);
microdosimetry->GetEntry(3);

-----
//how to make a TEntryList
TString filepath = "/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm.root";	

TFile *f = TFile::Open(filepath);
TTreeReader microdosimetryreader("microdosimetry", f);
TTreeReaderValue<Int_t> eventID(microdosimetryreader, "eventID");
TEntryList microdosimetry_entry_list = TEntryList();

while (microdosimetryreader.Next())
{
	if (*eventID < 1000)
	{
		microdosimetry_entry_list.Enter(microdosimetryreader.GetCurrentEntry());
	}
	else {break;}
}



TTreeReaderValue<Double_t> x(microdosimetryreader, "x");




TTree *microdosimetry;
f->GetObject("microdosimetry",microdosimetry);
//below takes a LONG time. There has to be a better way for variables that are already sorted.
microdosimetry->Draw("eventID","eventID<1000"); //draw only the first 1000 events

TString filepath = "/home/joseph/Documents/M1_local/trackfiles_1mil_livermore/co60_1cm.root";	
TFile *f = TFile::Open(filepath);
TTree *microdosimetry;
f->GetObject("microdosimetry",microdosimetry);

//Now, loop over the Tree cluster by cluster,
//and determine the entry and exit points for each TreeReader
Long64_t nentries = microdosimetry->GetEntries();
auto clusteriterator = microdosimetry->GetClusterIterator();

*/

