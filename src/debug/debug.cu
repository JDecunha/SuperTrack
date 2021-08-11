//Functions for debugging / verifying proper operation of the code

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