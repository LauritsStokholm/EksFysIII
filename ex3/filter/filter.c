// Read in digitiser data in list mode format and produce a TTree containing coincidence data

// Function to read in *.dat files with a graphical interface
void select_file(TString & filename, int opt){
    TGFileInfo file_info_;
    char *filetypes;
    if (opt==1){
        char *selection[] = {"List mode data", "*.txt", 0, 0};
        filetypes = selection;
    }
    else if (opt==2){
        char *selection[] = {"Calibration data", "*.cal", 0, 0};
        filetypes = selection;
    }
    else if (opt==3){
        char *selection[] = {"ROOT file", "*.root", 0, 0};
        filetypes = selection;
    }
    else if (opt==4){
        char *selection[] = {"Text file", "*.txt", 0, 0};
        filetypes = selection;
    }
    file_info_.fFileTypes = filetypes;
    file_info_.fIniDir = StrDup(".");
    
    new TGFileDialog(gClient->GetDefaultRoot(), gClient->GetDefaultRoot(), kFDOpen, &file_info_);
    if (file_info_.fFilename){
        cout << "'" << file_info_.fFilename << "' selected." << endl;
        filename = file_info_.fFilename;
    }
    else{
        cout << "Error: no file selected" << endl;
        exit(1);
    }
}

//Alternate readline function
int read_data(FILE *&infile, ULong64_t &t, Double_t &e, unsigned int &z, Double_t a, Double_t b){
    int test = fscanf(infile, "%llu %lf %d",&t,&e,&z);
    if (test == -1){ // EOF
        return 2;
    }
    else if (test != 3){ // Not 3 elements, probably header
        char buff[100];
        fgets(buff,100,infile); // Read the line into a temporary buffer instead
        return 1;
    }
    // Random number generator for smoothing calibrated energies
    TRandom *R = new TRandom(time(0));
    Double_t rndm = R->Rndm();
    e = (e+rndm-0.5)*a + b; // Calibrate E
    //cout << t << endl;
    return 0;
}

// Main function
void TTreeBuilder(){
    
    // Addresses of the two input files and one calibration file
    TString File0, File1, CalFile;
    
    cout << "Select input file for detector 0" << endl;
    select_file(File0,1);
    cout << "Select input file for detector 1" << endl;
    select_file(File1,1);
    cout << "Select calibration file" << endl;
    select_file(CalFile,2);
	
	// Read in the input files
	// This crashes on Windows if we don't convert to a char* first
	FILE *infile0, *infile1;
	char *filename0 = File0;
	char *filename1 = File1;
	
	infile0 = fopen( filename0, "r" );
	infile1 = fopen( filename1, "r" );

    // Address of an additional file to write raw ascii data to
    FILE *AsciiOut;
    TString AsciiFileName;
    cout << "Specify ASCII-format output file" << endl;
	select_file(AsciiFileName,4);
	char *filename3 = AsciiFileName;
    AsciiOut = fopen( filename3, "w");
    
    // Load in the paramaters from the calibration file
    gROOT->ProcessLine(".L " + CalFile);
    
    // Define a TTree to store our data and set up branches
    TTree *data = new TTree("data","data");

	TH1F * h0 = new TH1F("h0", "h0", 4000, 0., 16000.);
	TH1F * h1 = new TH1F("h1", "h1", 4000, 0., 16000.);
    
    Double_t E0, E1; // Energies
    ULong64_t T0, T1; // Time stamps
    Double_t TDiff; // Time difference for coincidence events
    unsigned int Z0, Z1; // Extras value
    
    data->Branch("E0",&E0,"E0/D"); data->Branch("E1",&E1,"E1/D");
    data->Branch("T0",&T0,"T0/l"); data->Branch("T1",&T1,"T1/l");
    data->Branch("TDiff",&TDiff,"TDiff/D");
    data->Branch("Z0",&Z0,"Z0/i"); data->Branch("Z1",&Z1,"Z1/i");
    
    // Bit to switch when a file ends
    bool EOF = false;
    
    // Count events
    int evts = 0;
    
    // Skip past headers and read in first event
    while (read_data(infile0,T0,E0,Z0,a0,b0) == 1);
    while (read_data(infile1,T1,E1,Z1,a1,b1) == 1);
    
    // Random number generator for smoothing
    TRandom *R = new TRandom(time(0));
    Double_t rndm;
    
    // Now start looping over data
    while(true){
        
        // At start of loop we should always have data
        rndm = R->Rndm();
        TDiff = Double_t(T0) - Double_t(T1) - TOff + rndm -0.5;
        
        // If TDiff<max, write and read new data
        while (abs(TDiff)>TWindow){
        
            if (T0<T1){
                if(read_data(infile0,T0,E0,Z0,a0,b0) == 2){
                    // EOF
                    EOF = true;
                    break;
                }
                rndm = R->Rndm();
                TDiff = Double_t(T0) - Double_t(T1) - TOff + rndm -0.5;
            }
            else if (T1<T0){
                if(read_data(infile1,T1,E1,Z1,a1,b1) == 2){
                    // EOF
                    EOF = true;
                    break;
                }
                rndm = R->Rndm();
                TDiff = Double_t(T0) - Double_t(T1) - TOff + rndm -0.5;
            }
        
        }
        
        if (EOF == true) break;
        
        // Coincidence data - write to file and read in new data from both files
        data->Fill();
        evts ++;
        if(evts%10000 == 0){
            cout << evts << " events processed\r" << flush;
            if (evts > (EvtLimit-10000) && EvtLimit>0){
                cout << "\n Event limit reached before end of file" << endl;
                break;
            }
        }

		h0->Fill(E0);
		h1->Fill(E1);
        
        // Also write to the ascii file
        fprintf(AsciiOut, "%llu %lf %llu %lf \n",T0, E0, T1, E1);
        
        if(read_data(infile0,T0,E0,Z0,a0,b0) == 2) break;
        if(read_data(infile1,T1,E1,Z1,a1,b1) == 2) break;
        
    }
    
    cout << "\nA total of " << evts << " coincidence events were found" << endl;
    
    // Two aliases to give time in seconds
    data->SetAlias("T0s","T0/1e8");
    data->SetAlias("T1s","T1/1e8");
    
    // Use another dialog to choose the save location for the final data
    cout << "Specify ROOT-format output file" << endl;
    TString OutFile;
    select_file(OutFile,3);
    TFile *tf = new TFile(OutFile,"RECREATE");
    data->Write();
	h0->Write();
	h1->Write();
    tf->Close();

}
