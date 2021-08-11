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

/*

Some various ROOT command line code for opening things and testing:


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

