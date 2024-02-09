#include "call_libraries.h" // call libraries from ROOT and C++
#include "uiclogo.h"	    // call UIC logo and initialization
#include "ntrkoff.h"        // get Ntrk offline

std::map<unsigned long long, int> runLumiEvtToEntryMap;
unsigned long long keyFromRunLumiEvent(UInt_t run, UInt_t lumi, ULong64_t event);

/*
Main skim pPb data and MC

Written by Dener Lemos (dener.lemos@cern.ch)

--> Arguments
input_file: text file with a list of root input files: Forest or Skims from jets
input_V0file: text files with a list of V0 files compatible with jets
ouputfile: just a counting number to run on Condor
*/
void V0Jet_pPbSkim_ZDC(TString input_file, TString input_V0file, TString ouputfile){

	float V0ptmin = 0.3;
	float V0etamin = 2.4;

	TString outputFileName;
	outputFileName = Form("%s",ouputfile.Data());

	clock_t sec_start, sec_end;
	sec_start = clock(); // start timing measurement

	TDatime* date = new TDatime();

	printwelcome(true); // welcome message
	print_start(); // start timing print

	// Read the input jet file(s)
	fstream inputfile;
	inputfile.open(Form("%s",input_file.Data()), ios::in);
	if(!inputfile.is_open()){cout << "List of JetTrack input files not founded!" << endl; return;}{cout << "List of JetTrack input files founded! --> " << input_file.Data() << endl;}

	std::vector<TString> file_name_vector;
	string file_chain;
	while(getline(inputfile, file_chain)){file_name_vector.push_back(Form("%s",file_chain.c_str()));}
	inputfile.close();

	TChain *heavyIonTree = new TChain("hiEvtAnalyzer/HiTree");
	TChain *skimTree = new TChain("skimanalysis/HltTree");
	TChain *trackTree = new TChain("ppTrack/trackTree");
	// add all the trees to the chain
	for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++){
		cout << "Adding file " << *listIterator << " to the chains" << endl;
		heavyIonTree->Add(*listIterator);
		skimTree->Add(*listIterator);
		trackTree->Add(*listIterator);
	}

	// Read the input V0 file(s)
	fstream inputfileV0;
	inputfileV0.open(Form("%s",input_V0file.Data()), ios::in);
	if(!inputfileV0.is_open()){cout << "List of V0 input files not founded!" << endl; return;}{cout << "List of V0 input files founded! --> " << input_V0file.Data() << endl;}

	// Make a chain and a vector of file names
	std::vector<TString> file_name_vectorV0;
	string file_chainV0;
	while(getline(inputfileV0, file_chainV0)){file_name_vectorV0.push_back(Form("%s",file_chainV0.c_str()));}
	inputfileV0.close();
	TChain *MainV0Tree = new TChain("K0SAnalysis/my_tree");
	TChain *K0sTree = new TChain("K0SAnalysis/my_treeK0s");
	TChain *LamTree = new TChain("K0SAnalysis/my_treeLam");
	TChain *CasTree = new TChain("K0SAnalysis/my_treeXi");
	TChain *OmeTree = new TChain("K0SAnalysis/my_treeOm");
	// add all the trees to the chain
	for (std::vector<TString>::iterator listIterator = file_name_vectorV0.begin(); listIterator != file_name_vectorV0.end(); listIterator++){
		cout << "Adding file " << *listIterator << " to the chains" << endl;
		MainV0Tree->Add(*listIterator);
		K0sTree->Add(*listIterator);
		LamTree->Add(*listIterator);
		CasTree->Add(*listIterator);
		OmeTree->Add(*listIterator);
	}

	// V0 Branches
	TBranch *V0_runBranch;	 // Branch for run
	TBranch *V0_evtBranch;	 // Branch for event
	TBranch *V0_lumiBranch;	 // Branch for lumi
	Int_t V0_run;			 // Run number
	Int_t V0_evt;			 // Event number
	Int_t V0_lumi;			 // Luminosity block

	//---------------------------  K0s  -------------------------
	//daughter 1 (pi+)
	TBranch *K0s_dxy1Branch; 
	TBranch *K0s_dz1Branch; 
	TBranch *K0s_chi21Branch; 
	TBranch *K0s_d1pxBranch; 
	TBranch *K0s_d1pyBranch; 
	TBranch *K0s_d1pzBranch; 
	TBranch *K0s_d1MBranch; 
	TBranch *K0s_d1NhitBranch;

	//daughter 2 (pi-)
	TBranch *K0s_dxy2Branch; 
	TBranch *K0s_dz2Branch; 
	TBranch *K0s_chi22Branch; 
	TBranch *K0s_d2pxBranch; 
	TBranch *K0s_d2pyBranch; 
	TBranch *K0s_d2pzBranch; 
	TBranch *K0s_d2MBranch; 
	TBranch *K0s_d2NhitBranch;

	//Mother (K0s)
	TBranch *K0s_3DaglBranch; 
	TBranch *K0s_3DdlBranch; 
	TBranch *K0s_ptBranch;
	TBranch *K0s_etaBranch; 
	TBranch *K0s_phiBranch; 
	TBranch *K0s_massBranch; 
	TBranch *K0s_dcaBranch; 
	TBranch *K0s_vtxBranch; 	

	// Leaves for K0s particles
	//daughter 1 (pi+)
	vector<double> *K0s_dxy1; 
	vector<double> *K0s_dz1; 
	vector<double> *K0s_chi21; 
	vector<double> *K0s_d1px; 
	vector<double> *K0s_d1py; 
	vector<double> *K0s_d1pz; 
	vector<double> *K0s_d1M; 
	vector<double> *K0s_d1Nhit;

	//daughter 2 (pi-)
	vector<double> *K0s_dxy2; 
	vector<double> *K0s_dz2; 
	vector<double> *K0s_chi22; 
	vector<double> *K0s_d2px; 
	vector<double> *K0s_d2py; 
	vector<double> *K0s_d2pz; 
	vector<double> *K0s_d2M; 
	vector<double> *K0s_d2Nhit;

	//Mother (K0s)
	vector<double> *K0s_3Dagl; 
	vector<double> *K0s_3Ddl; 
	vector<double> *K0s_pt; 
	vector<double> *K0s_eta; 
	vector<double> *K0s_phi; 
	vector<double> *K0s_mass; 
	vector<double> *K0s_dca; 
	vector<double> *K0s_vtx; 

	//---------------------------  Lambdas  -------------------------
	//daughter 1 (pi+)
	TBranch *Lam_dxy1Branch; 
	TBranch *Lam_dz1Branch; 
	TBranch *Lam_chi21Branch; 
	TBranch *Lam_d1pxBranch; 
	TBranch *Lam_d1pyBranch; 
	TBranch *Lam_d1pzBranch; 
	TBranch *Lam_d1MBranch; 
	TBranch *Lam_d1NhitBranch;

	//daughter 2 (pi-)
	TBranch *Lam_dxy2Branch; 
	TBranch *Lam_dz2Branch; 
	TBranch *Lam_chi22Branch; 
	TBranch *Lam_d2pxBranch; 
	TBranch *Lam_d2pyBranch; 
	TBranch *Lam_d2pzBranch; 
	TBranch *Lam_d2MBranch; 
	TBranch *Lam_d2NhitBranch;

	//Mother (Lambda)
	TBranch *Lam_3DaglBranch; 
	TBranch *Lam_3DdlBranch; 
	TBranch *Lam_ptBranch;
	TBranch *Lam_etaBranch; 
	TBranch *Lam_phiBranch; 
	TBranch *Lam_massBranch; 
	TBranch *Lam_dcaBranch; 
	TBranch *Lam_vtxBranch; 	
	TBranch *Lam_idBranch; 

	// Leaves for Lambda particles
	//daughter 1 (pi+)
	vector<double> *Lam_dxy1; 
	vector<double> *Lam_dz1; 
	vector<double> *Lam_chi21; 
	vector<double> *Lam_d1px; 
	vector<double> *Lam_d1py; 
	vector<double> *Lam_d1pz; 
	vector<double> *Lam_d1M; 
	vector<double> *Lam_d1Nhit;

	//daughter 2 (pi-)
	vector<double> *Lam_dxy2; 
	vector<double> *Lam_dz2; 
	vector<double> *Lam_chi22; 
	vector<double> *Lam_d2px; 
	vector<double> *Lam_d2py; 
	vector<double> *Lam_d2pz; 
	vector<double> *Lam_d2M; 
	vector<double> *Lam_d2Nhit;

	//Mother (Lambda)
	vector<double> *Lam_3Dagl; 
	vector<double> *Lam_3Ddl; 
	vector<double> *Lam_pt; 
	vector<double> *Lam_eta; 
	vector<double> *Lam_phi; 
	vector<double> *Lam_mass; 
	vector<double> *Lam_dca; 
	vector<double> *Lam_vtx; 
	vector<double> *Lam_id; 

	//---------------------------  Xi  -------------------------

	//daughter 1 (Lambda)
	TBranch *Xi_d1ptBranch;
	TBranch *Xi_d1etaBranch;
	TBranch *Xi_d1phiBranch;
	TBranch *Xi_d1massBranch;

	//Lambda daughters
	//proton
	TBranch *Xi_chi21_1Branch;
	TBranch *Xi_d1pt_1Branch;
	TBranch *Xi_d1eta_1Branch;
	TBranch *Xi_d1phi_1Branch;
	TBranch *Xi_d1mass_1Branch;
	TBranch *Xi_d1Nhit_1Branch;

	//pion
	TBranch *Xi_chi21_2Branch;
	TBranch *Xi_d1pt_2Branch;
	TBranch *Xi_d1eta_2Branch;
	TBranch *Xi_d1phi_2Branch;
	TBranch *Xi_d1mass_2Branch;
	TBranch *Xi_d1Nhit_2Branch;

	//daughter 2 (pion)
	TBranch *Xi_chi22Branch;
	TBranch *Xi_d2ptBranch;
	TBranch *Xi_d2etaBranch;
	TBranch *Xi_d2phiBranch;
	TBranch *Xi_d2massBranch;
	TBranch *Xi_d2NhitBranch;

	//Mother (Xi)
	TBranch *Xi_cas3DIpSigValueBranch;
	TBranch *Xi_casPi3DIpSigValueBranch;
	TBranch *Xi_VTrkPi3DIpSigValueBranch;
	TBranch *Xi_VTrkP3DIpSigValueBranch;
	TBranch *Xi_casFlightSigValueBranch;
	TBranch *Xi_distanceSigValueBranch;
	TBranch *Xi_ptBranch;
	TBranch *Xi_etaBranch;
	TBranch *Xi_phiBranch;
	TBranch *Xi_massBranch;
	TBranch *Xi_idBranch;

	//daughter 1 (Lambda)
	vector<double> *Xi_d1pt;
	vector<double> *Xi_d1eta;
	vector<double> *Xi_d1phi;
	vector<double> *Xi_d1mass;

	//Lambda daughters
	//proton
	vector<double> *Xi_chi21_1;
	vector<double> *Xi_d1pt_1;
	vector<double> *Xi_d1eta_1;
	vector<double> *Xi_d1phi_1;
	vector<double> *Xi_d1mass_1;
	vector<double> *Xi_d1Nhit_1;

	//pion
	vector<double> *Xi_chi21_2;
	vector<double> *Xi_d1pt_2;
	vector<double> *Xi_d1eta_2;
	vector<double> *Xi_d1phi_2;
	vector<double> *Xi_d1mass_2;
	vector<double> *Xi_d1Nhit_2;

	//daughter 2 (pion)
	vector<double> *Xi_chi22;
	vector<double> *Xi_d2pt;
	vector<double> *Xi_d2eta;
	vector<double> *Xi_d2phi;
	vector<double> *Xi_d2mass;
	vector<double> *Xi_d2Nhit;

	//Mother (Xi)
	vector<double> *Xi_cas3DIpSigValue;
	vector<double> *Xi_casPi3DIpSigValue;
	vector<double> *Xi_VTrkPi3DIpSigValue;
	vector<double> *Xi_VTrkP3DIpSigValue;
	vector<double> *Xi_casFlightSigValue;
	vector<double> *Xi_distanceSigValue;
	vector<double> *Xi_pt;
	vector<double> *Xi_eta;
	vector<double> *Xi_phi;
	vector<double> *Xi_mass;
	vector<double> *Xi_id;

	//---------------------------  Omega  -------------------------

	//daughter 1 (Lambda)
	TBranch *Om_d1ptBranch;
	TBranch *Om_d1etaBranch;
	TBranch *Om_d1phiBranch;
	TBranch *Om_d1massBranch;

	//Lambda daughters
	//proton
	TBranch *Om_chi21_1Branch;
	TBranch *Om_d1pt_1Branch;
	TBranch *Om_d1eta_1Branch;
	TBranch *Om_d1phi_1Branch;
	TBranch *Om_d1mass_1Branch;
	TBranch *Om_d1Nhit_1Branch;

	//pion
	TBranch *Om_chi21_2Branch;
	TBranch *Om_d1pt_2Branch;
	TBranch *Om_d1eta_2Branch;
	TBranch *Om_d1phi_2Branch;
	TBranch *Om_d1mass_2Branch;
	TBranch *Om_d1Nhit_2Branch;

	//daughter 2 (pion)
	TBranch *Om_chi22Branch;
	TBranch *Om_d2ptBranch;
	TBranch *Om_d2etaBranch;
	TBranch *Om_d2phiBranch;
	TBranch *Om_d2massBranch;
	TBranch *Om_d2NhitBranch;

	//Mother (Om)
	TBranch *Om_cas3DIpSigValueBranch;
	TBranch *Om_casPi3DIpSigValueBranch;
	TBranch *Om_VTrkPi3DIpSigValueBranch;
	TBranch *Om_VTrkP3DIpSigValueBranch;
	TBranch *Om_casFlightSigValueBranch;
	TBranch *Om_distanceSigValueBranch;
	TBranch *Om_ptBranch;
	TBranch *Om_etaBranch;
	TBranch *Om_phiBranch;
	TBranch *Om_massBranch;
	TBranch *Om_idBranch;

	//daughter 1 (Lambda)
	vector<double> *Om_d1pt;
	vector<double> *Om_d1eta;
	vector<double> *Om_d1phi;
	vector<double> *Om_d1mass;

	//Lambda daughters
	//proton
	vector<double> *Om_chi21_1;
	vector<double> *Om_d1pt_1;
	vector<double> *Om_d1eta_1;
	vector<double> *Om_d1phi_1;
	vector<double> *Om_d1mass_1;
	vector<double> *Om_d1Nhit_1;

	//pion
	vector<double> *Om_chi21_2;
	vector<double> *Om_d1pt_2;
	vector<double> *Om_d1eta_2;
	vector<double> *Om_d1phi_2;
	vector<double> *Om_d1mass_2;
	vector<double> *Om_d1Nhit_2;

	//daughter 2 (pion)
	vector<double> *Om_chi22;
	vector<double> *Om_d2pt;
	vector<double> *Om_d2eta;
	vector<double> *Om_d2phi;
	vector<double> *Om_d2mass;
	vector<double> *Om_d2Nhit;

	//Mother (Om)
	vector<double> *Om_cas3DIpSigValue;
	vector<double> *Om_casPi3DIpSigValue;
	vector<double> *Om_VTrkPi3DIpSigValue;
	vector<double> *Om_VTrkP3DIpSigValue;
	vector<double> *Om_casFlightSigValue;
	vector<double> *Om_distanceSigValue;
	vector<double> *Om_pt;
	vector<double> *Om_eta;
	vector<double> *Om_phi;
	vector<double> *Om_mass;
	vector<double> *Om_id;

	//Main tree (event info)
	MainV0Tree->SetBranchStatus("*",0);
	MainV0Tree->SetBranchStatus("runNumber",1);
	MainV0Tree->SetBranchAddress("runNumber",&V0_run,&V0_runBranch);
	MainV0Tree->SetBranchStatus("evNumber",1);
	MainV0Tree->SetBranchAddress("evNumber",&V0_evt,&V0_evtBranch);
	MainV0Tree->SetBranchStatus("LumiSection",1);
	MainV0Tree->SetBranchAddress("LumiSection",&V0_lumi,&V0_lumiBranch);

	//K0s tree
	K0sTree->SetBranchStatus("*",0);

	K0sTree->SetBranchStatus("K0s_dxy1",1);
	K0sTree->SetBranchAddress("K0s_dxy1",&K0s_dxy1,&K0s_dxy1Branch);
	K0sTree->SetBranchStatus("K0s_dz1",1);
	K0sTree->SetBranchAddress("K0s_dz1",&K0s_dz1,&K0s_dz1Branch);
	K0sTree->SetBranchStatus("K0s_chi21",1);
	K0sTree->SetBranchAddress("K0s_chi21",&K0s_chi21,&K0s_chi21Branch);
	K0sTree->SetBranchStatus("K0s_d1px",1);
	K0sTree->SetBranchAddress("K0s_d1px",&K0s_d1px,&K0s_d1pxBranch);
	K0sTree->SetBranchStatus("K0s_d1py",1);
	K0sTree->SetBranchAddress("K0s_d1py",&K0s_d1py,&K0s_d1pyBranch);
	K0sTree->SetBranchStatus("K0s_d1pz",1);
	K0sTree->SetBranchAddress("K0s_d1pz",&K0s_d1pz,&K0s_d1pzBranch);
	K0sTree->SetBranchStatus("K0s_d1M",1);
	K0sTree->SetBranchAddress("K0s_d1M",&K0s_d1M,&K0s_d1MBranch);
	K0sTree->SetBranchStatus("K0s_d1Nhit",1);
	K0sTree->SetBranchAddress("K0s_d1Nhit",&K0s_d1Nhit,&K0s_d1NhitBranch);

	K0sTree->SetBranchStatus("K0s_dxy2",1);
	K0sTree->SetBranchAddress("K0s_dxy2",&K0s_dxy2,&K0s_dxy2Branch);
	K0sTree->SetBranchStatus("K0s_dz2",1);
	K0sTree->SetBranchAddress("K0s_dz2",&K0s_dz2,&K0s_dz2Branch);
	K0sTree->SetBranchStatus("K0s_chi22",1);
	K0sTree->SetBranchAddress("K0s_chi22",&K0s_chi22,&K0s_chi22Branch);
	K0sTree->SetBranchStatus("K0s_d2px",1);
	K0sTree->SetBranchAddress("K0s_d2px",&K0s_d2px,&K0s_d2pxBranch);
	K0sTree->SetBranchStatus("K0s_d2py",1);
	K0sTree->SetBranchAddress("K0s_d2py",&K0s_d2py,&K0s_d2pyBranch);
	K0sTree->SetBranchStatus("K0s_d2pz",1);
	K0sTree->SetBranchAddress("K0s_d2pz",&K0s_d2pz,&K0s_d2pzBranch);
	K0sTree->SetBranchStatus("K0s_d2M",1);
	K0sTree->SetBranchAddress("K0s_d2M",&K0s_d2M,&K0s_d2MBranch);
	K0sTree->SetBranchStatus("K0s_d2Nhit",1);
	K0sTree->SetBranchAddress("K0s_d2Nhit",&K0s_d2Nhit,&K0s_d2NhitBranch);

	K0sTree->SetBranchStatus("K0s_3Dagl",1);
	K0sTree->SetBranchAddress("K0s_3Dagl",&K0s_3Dagl,&K0s_3DaglBranch);	
	K0sTree->SetBranchStatus("K0s_3Ddl",1);
	K0sTree->SetBranchAddress("K0s_3Ddl",&K0s_3Ddl,&K0s_3DdlBranch);	
	K0sTree->SetBranchStatus("K0s_pt",1);
	K0sTree->SetBranchAddress("K0s_pt",&K0s_pt,&K0s_ptBranch);	
	K0sTree->SetBranchStatus("K0s_eta",1);
	K0sTree->SetBranchAddress("K0s_eta",&K0s_eta,&K0s_etaBranch);	
	K0sTree->SetBranchStatus("K0s_phi",1);
	K0sTree->SetBranchAddress("K0s_phi",&K0s_phi,&K0s_phiBranch);	
	K0sTree->SetBranchStatus("K0s_mass",1);
	K0sTree->SetBranchAddress("K0s_mass",&K0s_mass,&K0s_massBranch);	
	K0sTree->SetBranchStatus("K0s_dca",1);
	K0sTree->SetBranchAddress("K0s_dca",&K0s_dca,&K0s_dcaBranch);	
	K0sTree->SetBranchStatus("K0s_vtx",1);
	K0sTree->SetBranchAddress("K0s_vtx",&K0s_vtx,&K0s_vtxBranch);	

	//Lam tree

	LamTree->SetBranchStatus("*",0);

	LamTree->SetBranchStatus("Lam_dxy1",1);
	LamTree->SetBranchAddress("Lam_dxy1",&Lam_dxy1,&Lam_dxy1Branch);
	LamTree->SetBranchStatus("Lam_dz1",1);
	LamTree->SetBranchAddress("Lam_dz1",&Lam_dz1,&Lam_dz1Branch);
	LamTree->SetBranchStatus("Lam_chi21",1);
	LamTree->SetBranchAddress("Lam_chi21",&Lam_chi21,&Lam_chi21Branch);
	LamTree->SetBranchStatus("Lam_d1px",1);
	LamTree->SetBranchAddress("Lam_d1px",&Lam_d1px,&Lam_d1pxBranch);
	LamTree->SetBranchStatus("Lam_d1py",1);
	LamTree->SetBranchAddress("Lam_d1py",&Lam_d1py,&Lam_d1pyBranch);
	LamTree->SetBranchStatus("Lam_d1pz",1);
	LamTree->SetBranchAddress("Lam_d1pz",&Lam_d1pz,&Lam_d1pzBranch);
	LamTree->SetBranchStatus("Lam_d1M",1);
	LamTree->SetBranchAddress("Lam_d1M",&Lam_d1M,&Lam_d1MBranch);
	LamTree->SetBranchStatus("Lam_d1Nhit",1);
	LamTree->SetBranchAddress("Lam_d1Nhit",&Lam_d1Nhit,&Lam_d1NhitBranch);

	LamTree->SetBranchStatus("Lam_dxy2",1);
	LamTree->SetBranchAddress("Lam_dxy2",&Lam_dxy2,&Lam_dxy2Branch);
	LamTree->SetBranchStatus("Lam_dz2",1);
	LamTree->SetBranchAddress("Lam_dz2",&Lam_dz2,&Lam_dz2Branch);
	LamTree->SetBranchStatus("Lam_chi22",1);
	LamTree->SetBranchAddress("Lam_chi22",&Lam_chi22,&Lam_chi22Branch);
	LamTree->SetBranchStatus("Lam_d2px",1);
	LamTree->SetBranchAddress("Lam_d2px",&Lam_d2px,&Lam_d2pxBranch);
	LamTree->SetBranchStatus("Lam_d2py",1);
	LamTree->SetBranchAddress("Lam_d2py",&Lam_d2py,&Lam_d2pyBranch);
	LamTree->SetBranchStatus("Lam_d2pz",1);
	LamTree->SetBranchAddress("Lam_d2pz",&Lam_d2pz,&Lam_d2pzBranch);
	LamTree->SetBranchStatus("Lam_d2M",1);
	LamTree->SetBranchAddress("Lam_d2M",&Lam_d2M,&Lam_d2MBranch);
	LamTree->SetBranchStatus("Lam_d2Nhit",1);
	LamTree->SetBranchAddress("Lam_d2Nhit",&Lam_d2Nhit,&Lam_d2NhitBranch);

	LamTree->SetBranchStatus("Lam_3Dagl",1);
	LamTree->SetBranchAddress("Lam_3Dagl",&Lam_3Dagl,&Lam_3DaglBranch);	
	LamTree->SetBranchStatus("Lam_3Ddl",1);
	LamTree->SetBranchAddress("Lam_3Ddl",&Lam_3Ddl,&Lam_3DdlBranch);	
	LamTree->SetBranchStatus("Lam_pt",1);
	LamTree->SetBranchAddress("Lam_pt",&Lam_pt,&Lam_ptBranch);	
	LamTree->SetBranchStatus("Lam_eta",1);
	LamTree->SetBranchAddress("Lam_eta",&Lam_eta,&Lam_etaBranch);	
	LamTree->SetBranchStatus("Lam_phi",1);
	LamTree->SetBranchAddress("Lam_phi",&Lam_phi,&Lam_phiBranch);	
	LamTree->SetBranchStatus("Lam_mass",1);
	LamTree->SetBranchAddress("Lam_mass",&Lam_mass,&Lam_massBranch);	
	LamTree->SetBranchStatus("Lam_dca",1);
	LamTree->SetBranchAddress("Lam_dca",&Lam_dca,&Lam_dcaBranch);	
	LamTree->SetBranchStatus("Lam_vtx",1);
	LamTree->SetBranchAddress("Lam_vtx",&Lam_vtx,&Lam_vtxBranch);	
	LamTree->SetBranchStatus("Lam_id",1);
	LamTree->SetBranchAddress("Lam_id",&Lam_id,&Lam_idBranch);	

	CasTree->SetBranchStatus("Xi_d1pt",1);
	CasTree->SetBranchAddress("Xi_d1pt",&Xi_d1pt,&Xi_d1ptBranch);
	CasTree->SetBranchStatus("Xi_d1eta",1);
	CasTree->SetBranchAddress("Xi_d1eta",&Xi_d1eta,&Xi_d1etaBranch);
	CasTree->SetBranchStatus("Xi_d1phi",1);
	CasTree->SetBranchAddress("Xi_d1phi",&Xi_d1phi,&Xi_d1phiBranch);
	CasTree->SetBranchStatus("Xi_d1mass",1);
	CasTree->SetBranchAddress("Xi_d1mass",&Xi_d1mass,&Xi_d1massBranch);

	CasTree->SetBranchStatus("Xi_chi21_1",1);
	CasTree->SetBranchAddress("Xi_chi21_1",&Xi_chi21_1,&Xi_chi21_1Branch);
	CasTree->SetBranchStatus("Xi_d1pt_1",1);
	CasTree->SetBranchAddress("Xi_d1pt_1",&Xi_d1pt_1,&Xi_d1pt_1Branch);
	CasTree->SetBranchStatus("Xi_d1eta_1",1);
	CasTree->SetBranchAddress("Xi_d1eta_1",&Xi_d1eta_1,&Xi_d1eta_1Branch);
	CasTree->SetBranchStatus("Xi_d1phi_1",1);
	CasTree->SetBranchAddress("Xi_d1phi_1",&Xi_d1phi_1,&Xi_d1phi_1Branch);
	CasTree->SetBranchStatus("Xi_d1mass_1",1);
	CasTree->SetBranchAddress("Xi_d1mass_1",&Xi_d1mass_1,&Xi_d1mass_1Branch);
	CasTree->SetBranchStatus("Xi_d1Nhit_1",1);
	CasTree->SetBranchAddress("Xi_d1Nhit_1",&Xi_d1Nhit_1,&Xi_d1Nhit_1Branch);

	CasTree->SetBranchStatus("Xi_chi21_2",1);
	CasTree->SetBranchAddress("Xi_chi21_2",&Xi_chi21_2,&Xi_chi21_2Branch);
	CasTree->SetBranchStatus("Xi_d1pt_2",1);
	CasTree->SetBranchAddress("Xi_d1pt_2",&Xi_d1pt_2,&Xi_d1pt_2Branch);
	CasTree->SetBranchStatus("Xi_d1eta_2",1);
	CasTree->SetBranchAddress("Xi_d1eta_2",&Xi_d1eta_2,&Xi_d1eta_2Branch);
	CasTree->SetBranchStatus("Xi_d1phi_2",1);
	CasTree->SetBranchAddress("Xi_d1phi_2",&Xi_d1phi_2,&Xi_d1phi_2Branch);
	CasTree->SetBranchStatus("Xi_d1mass_2",1);
	CasTree->SetBranchAddress("Xi_d1mass_2",&Xi_d1mass_2,&Xi_d1mass_2Branch);
	CasTree->SetBranchStatus("Xi_d1Nhit_2",1);
	CasTree->SetBranchAddress("Xi_d1Nhit_2",&Xi_d1Nhit_2,&Xi_d1Nhit_2Branch);

	CasTree->SetBranchStatus("Xi_chi22",1);
	CasTree->SetBranchAddress("Xi_chi22",&Xi_chi22,&Xi_chi22Branch);
	CasTree->SetBranchStatus("Xi_d2pt",1);
	CasTree->SetBranchAddress("Xi_d2pt",&Xi_d2pt,&Xi_d2ptBranch);
	CasTree->SetBranchStatus("Xi_d2eta",1);
	CasTree->SetBranchAddress("Xi_d2eta",&Xi_d2eta,&Xi_d2etaBranch);
	CasTree->SetBranchStatus("Xi_d2phi",1);
	CasTree->SetBranchAddress("Xi_d2phi",&Xi_d2phi,&Xi_d2phiBranch);
	CasTree->SetBranchStatus("Xi_d2mass",1);
	CasTree->SetBranchAddress("Xi_d2mass",&Xi_d2mass,&Xi_d2massBranch);
	CasTree->SetBranchStatus("Xi_d2Nhit",1);
	CasTree->SetBranchAddress("Xi_d2Nhit",&Xi_d2Nhit,&Xi_d2NhitBranch);	

	CasTree->SetBranchStatus("Xi_cas3DIpSigValue",1);
	CasTree->SetBranchAddress("Xi_cas3DIpSigValue",&Xi_cas3DIpSigValue,&Xi_cas3DIpSigValueBranch);	
	CasTree->SetBranchStatus("Xi_casPi3DIpSigValue",1);
	CasTree->SetBranchAddress("Xi_casPi3DIpSigValue",&Xi_casPi3DIpSigValue,&Xi_casPi3DIpSigValueBranch);	
	CasTree->SetBranchStatus("Xi_VTrkPi3DIpSigValue",1);
	CasTree->SetBranchAddress("Xi_VTrkPi3DIpSigValue",&Xi_VTrkPi3DIpSigValue,&Xi_VTrkPi3DIpSigValueBranch);	
	CasTree->SetBranchStatus("Xi_VTrkP3DIpSigValue",1);
	CasTree->SetBranchAddress("Xi_VTrkP3DIpSigValue",&Xi_VTrkP3DIpSigValue,&Xi_VTrkP3DIpSigValueBranch);	
	CasTree->SetBranchStatus("Xi_casFlightSigValue",1);
	CasTree->SetBranchAddress("Xi_casFlightSigValue",&Xi_casFlightSigValue,&Xi_casFlightSigValueBranch);
	CasTree->SetBranchStatus("Xi_distanceSigValue",1);
	CasTree->SetBranchAddress("Xi_distanceSigValue",&Xi_distanceSigValue,&Xi_distanceSigValueBranch);
	CasTree->SetBranchStatus("Xi_pt",1);
	CasTree->SetBranchAddress("Xi_pt",&Xi_pt,&Xi_ptBranch);	
	CasTree->SetBranchStatus("Xi_eta",1);
	CasTree->SetBranchAddress("Xi_eta",&Xi_eta,&Xi_etaBranch);	
	CasTree->SetBranchStatus("Xi_phi",1);
	CasTree->SetBranchAddress("Xi_phi",&Xi_phi,&Xi_phiBranch);	
	CasTree->SetBranchStatus("Xi_mass",1);
	CasTree->SetBranchAddress("Xi_mass",&Xi_mass,&Xi_massBranch);	
	CasTree->SetBranchStatus("Xi_id",1);
	CasTree->SetBranchAddress("Xi_id",&Xi_id,&Xi_idBranch);

	OmeTree->SetBranchStatus("Om_d1pt",1);
	OmeTree->SetBranchAddress("Om_d1pt",&Om_d1pt,&Om_d1ptBranch);
	OmeTree->SetBranchStatus("Om_d1eta",1);
	OmeTree->SetBranchAddress("Om_d1eta",&Om_d1eta,&Om_d1etaBranch);
	OmeTree->SetBranchStatus("Om_d1phi",1);
	OmeTree->SetBranchAddress("Om_d1phi",&Om_d1phi,&Om_d1phiBranch);
	OmeTree->SetBranchStatus("Om_d1mass",1);
	OmeTree->SetBranchAddress("Om_d1mass",&Om_d1mass,&Om_d1massBranch);

	OmeTree->SetBranchStatus("Om_chi21_1",1);
	OmeTree->SetBranchAddress("Om_chi21_1",&Om_chi21_1,&Om_chi21_1Branch);
	OmeTree->SetBranchStatus("Om_d1pt_1",1);
	OmeTree->SetBranchAddress("Om_d1pt_1",&Om_d1pt_1,&Om_d1pt_1Branch);
	OmeTree->SetBranchStatus("Om_d1eta_1",1);
	OmeTree->SetBranchAddress("Om_d1eta_1",&Om_d1eta_1,&Om_d1eta_1Branch);
	OmeTree->SetBranchStatus("Om_d1phi_1",1);
	OmeTree->SetBranchAddress("Om_d1phi_1",&Om_d1phi_1,&Om_d1phi_1Branch);
	OmeTree->SetBranchStatus("Om_d1mass_1",1);
	OmeTree->SetBranchAddress("Om_d1mass_1",&Om_d1mass_1,&Om_d1mass_1Branch);
	OmeTree->SetBranchStatus("Om_d1Nhit_1",1);
	OmeTree->SetBranchAddress("Om_d1Nhit_1",&Om_d1Nhit_1,&Om_d1Nhit_1Branch);

	OmeTree->SetBranchStatus("Om_chi21_2",1);
	OmeTree->SetBranchAddress("Om_chi21_2",&Om_chi21_2,&Om_chi21_2Branch);
	OmeTree->SetBranchStatus("Om_d1pt_2",1);
	OmeTree->SetBranchAddress("Om_d1pt_2",&Om_d1pt_2,&Om_d1pt_2Branch);
	OmeTree->SetBranchStatus("Om_d1eta_2",1);
	OmeTree->SetBranchAddress("Om_d1eta_2",&Om_d1eta_2,&Om_d1eta_2Branch);
	OmeTree->SetBranchStatus("Om_d1phi_2",1);
	OmeTree->SetBranchAddress("Om_d1phi_2",&Om_d1phi_2,&Om_d1phi_2Branch);
	OmeTree->SetBranchStatus("Om_d1mass_2",1);
	OmeTree->SetBranchAddress("Om_d1mass_2",&Om_d1mass_2,&Om_d1mass_2Branch);
	OmeTree->SetBranchStatus("Om_d1Nhit_2",1);
	OmeTree->SetBranchAddress("Om_d1Nhit_2",&Om_d1Nhit_2,&Om_d1Nhit_2Branch);

	OmeTree->SetBranchStatus("Om_chi22",1);
	OmeTree->SetBranchAddress("Om_chi22",&Om_chi22,&Om_chi22Branch);
	OmeTree->SetBranchStatus("Om_d2pt",1);
	OmeTree->SetBranchAddress("Om_d2pt",&Om_d2pt,&Om_d2ptBranch);
	OmeTree->SetBranchStatus("Om_d2eta",1);
	OmeTree->SetBranchAddress("Om_d2eta",&Om_d2eta,&Om_d2etaBranch);
	OmeTree->SetBranchStatus("Om_d2phi",1);
	OmeTree->SetBranchAddress("Om_d2phi",&Om_d2phi,&Om_d2phiBranch);
	OmeTree->SetBranchStatus("Om_d2mass",1);
	OmeTree->SetBranchAddress("Om_d2mass",&Om_d2mass,&Om_d2massBranch);
	OmeTree->SetBranchStatus("Om_d2Nhit",1);
	OmeTree->SetBranchAddress("Om_d2Nhit",&Om_d2Nhit,&Om_d2NhitBranch);	

	OmeTree->SetBranchStatus("Om_cas3DIpSigValue",1);
	OmeTree->SetBranchAddress("Om_cas3DIpSigValue",&Om_cas3DIpSigValue,&Om_cas3DIpSigValueBranch);	
	OmeTree->SetBranchStatus("Om_casPi3DIpSigValue",1);
	OmeTree->SetBranchAddress("Om_casPi3DIpSigValue",&Om_casPi3DIpSigValue,&Om_casPi3DIpSigValueBranch);	
	OmeTree->SetBranchStatus("Om_VTrkPi3DIpSigValue",1);
	OmeTree->SetBranchAddress("Om_VTrkPi3DIpSigValue",&Om_VTrkPi3DIpSigValue,&Om_VTrkPi3DIpSigValueBranch);	
	OmeTree->SetBranchStatus("Om_VTrkP3DIpSigValue",1);
	OmeTree->SetBranchAddress("Om_VTrkP3DIpSigValue",&Om_VTrkP3DIpSigValue,&Om_VTrkP3DIpSigValueBranch);	
	OmeTree->SetBranchStatus("Om_casFlightSigValue",1);
	OmeTree->SetBranchAddress("Om_casFlightSigValue",&Om_casFlightSigValue,&Om_casFlightSigValueBranch);
	OmeTree->SetBranchStatus("Om_distanceSigValue",1);
	OmeTree->SetBranchAddress("Om_distanceSigValue",&Om_distanceSigValue,&Om_distanceSigValueBranch);
	OmeTree->SetBranchStatus("Om_pt",1);
	OmeTree->SetBranchAddress("Om_pt",&Om_pt,&Om_ptBranch);	
	OmeTree->SetBranchStatus("Om_eta",1);
	OmeTree->SetBranchAddress("Om_eta",&Om_eta,&Om_etaBranch);	
	OmeTree->SetBranchStatus("Om_phi",1);
	OmeTree->SetBranchAddress("Om_phi",&Om_phi,&Om_phiBranch);	
	OmeTree->SetBranchStatus("Om_mass",1);
	OmeTree->SetBranchAddress("Om_mass",&Om_mass,&Om_massBranch);	
	OmeTree->SetBranchStatus("Om_id",1);
	OmeTree->SetBranchAddress("Om_id",&Om_id,&Om_idBranch);

	// Bellow here is for jets
	// Branches for heavy ion tree
	TBranch *runBranch;						 	// Branch for run
	TBranch *eventBranch;					 	// Branch for event
	TBranch *lumiBranch;						// Branch for lumi
	TBranch *hiVzBranch;						// Branch for vertex z-position
	TBranch *hiHFplusBranch;					// Branch for HF+ energy deposity
	TBranch *hiHFminusBranch;					// Branch for HF- energy deposity
	TBranch *hiZDCplusBranch;					// Branch for ZDC+ energy deposity
	TBranch *hiZDCminusBranch;					// Branch for ZDC- energy deposity
	TBranch *hi_FRGBranch;						// Branch for Forward Rapidity Gap

	// Leaves for heavy ion tree
	UInt_t run;					 // Run number
	ULong64_t event;			 // Event number
	UInt_t lumi;				 // Luminosity block
	Float_t vertexZ;			 // Vertex z-position
	Float_t hiHFplus;			 // transverse energy sum of HF+ tower;
	Float_t hiHFminus;			 // transverse energy sum of HF- tower;
	Float_t hiZDCplus;			 // transverse energy sum of HF+ tower;
	Float_t hiZDCminus;			 // transverse energy sum of HF- tower;
	Float_t hi_FRG;			 	 // Forward Rapidity Gap;
	Int_t Ntrkoff;			 	 // Offline multiplicity
	Float_t NtrkCorr;		 	 // Corrected multiplicity
	Float_t NtrkCorrTight;	 	 // Tight multiplicity
	Float_t NtrkCorrLoose;	 	 // Loose multiplicity

	// Branches for skim tree
	TBranch *primaryVertexBranch;						// Branch for primary vertex filter bit
	TBranch *beamScrapingBranch;				// Branch for beam scraping filter bit
	TBranch *hfCoincidenceBranch;				// Branch for energy recorded one HF tower above threshold on each side
	TBranch *pVertexFilterCutdz1p0Branch;		// Branch for PU Filter default
	TBranch *pVertexFilterCutGplusBranch;		// Branch for PU Filter GPlus
	TBranch *pVertexFilterCutVtx1Branch;		// Branch for PU Filter 1 vertex only

	// Leaves for the skim tree
	Int_t primaryVertexFilterBit;				// Filter bit for primary vertex
	Int_t beamScrapingFilterBit;				// Filter bit for beam scraping
	Int_t hfCoincidenceFilterBit;				// Filter bit or energy recorded one HF tower above threshold on each side
	Int_t pVertexFilterCutdz1p0Bit;				// Filter bit for PU Filter
	Int_t pVertexFilterCutGplusBit;				// Filter bit for PU Filter
	Int_t pVertexFilterCutVtx1Bit;					// Filter bit for PU Filter

	// Branches for track tree
	TBranch *nTracksBranch;									// Branch for number of tracks
	TBranch *trackPtBranch;									// Branch for track pT
	TBranch *trackPtErrorBranch;							// Branch for track pT error
	TBranch *trackPhiBranch;								// Branch for track phi
	TBranch *trackEtaBranch;								// Branch for track eta
	TBranch *trackHighPurityBranch;							// Branch for high purity of the track
	TBranch *trackVertexDistanceZBranch;			 		// Branch for track distance from primary vertex in z-direction
	TBranch *trackVertexDistanceZErrorBranch;				// Branch for error for track distance from primary vertex in z-direction
	TBranch *trackVertexDistanceXYBranch;					// Branch for track distance from primary vertex in xy-direction
	TBranch *trackVertexDistanceXYErrorBranch; 				// Branch for error for track distance from primary vertex in xy-direction
	TBranch *trackChargeBranch;								// Branch for track charge
	
	// Leaves for the track tree
	const Int_t nMaxTrack = 2000;
	Int_t nTracks;														// Number of tracks
	Float_t trackPtArray[nMaxTrack] = {0};								// Array for track pT
	Float_t trackPtErrorArray[nMaxTrack] = {0};							// Array for track pT errors
	Float_t trackPhiArray[nMaxTrack] = {0};								// Array for track phis
	Float_t trackEtaArray[nMaxTrack] = {0};								// Array for track etas
	Bool_t trackHighPurityArray[nMaxTrack] = {0};						// Array for the high purity of tracks
	Float_t trackVertexDistanceZArray[nMaxTrack] = {0};			 		// Array for track distance from primary vertex in z-direction
	Float_t trackVertexDistanceZErrorArray[nMaxTrack] = {0};			// Array for error for track distance from primary vertex in z-direction
	Float_t trackVertexDistanceXYArray[nMaxTrack] = {0};				// Array for track distance from primary vertex in xy-direction
	Float_t trackVertexDistanceXYErrorArray[nMaxTrack] = {0}; 			// Array for error for track distance from primary vertex in xy-direction
	Int_t trackChargeArray[nMaxTrack] = {0}; 										// Array for track charge


	// ========================================== //
	// Read all the branches from the input trees //
	// ========================================== //
	
	// Connect the branches of the heavy ion tree
	heavyIonTree->SetBranchStatus("*",0); // remove all branchs to read it fast
	heavyIonTree->SetBranchStatus("run",1);
	heavyIonTree->SetBranchAddress("run",&run,&runBranch);
	heavyIonTree->SetBranchStatus("evt",1);
	heavyIonTree->SetBranchAddress("evt",&event,&eventBranch);
	heavyIonTree->SetBranchStatus("lumi",1);
	heavyIonTree->SetBranchAddress("lumi",&lumi,&lumiBranch);
	heavyIonTree->SetBranchStatus("vz",1);
	heavyIonTree->SetBranchAddress("vz",&vertexZ,&hiVzBranch);
	heavyIonTree->SetBranchStatus("hiHFplus",1);
	heavyIonTree->SetBranchAddress("hiHFplus",&hiHFplus,&hiHFplusBranch);
	heavyIonTree->SetBranchStatus("hiHFminus",1);
	heavyIonTree->SetBranchAddress("hiHFminus",&hiHFminus,&hiHFminusBranch);
	heavyIonTree->SetBranchStatus("hiZDCplus",1);
	heavyIonTree->SetBranchAddress("hiZDCplus",&hiZDCplus,&hiZDCplusBranch);
	heavyIonTree->SetBranchStatus("hiZDCminus",1);
	heavyIonTree->SetBranchAddress("hiZDCminus",&hiZDCminus,&hiZDCminusBranch);
	heavyIonTree->SetBranchStatus("hi_FRG",1);
	heavyIonTree->SetBranchAddress("hi_FRG",&hi_FRG,&hi_FRGBranch);
		
	// Connect the branches to the skim tree
	skimTree->SetBranchStatus("*",0);
	skimTree->SetBranchStatus("pPAprimaryVertexFilter",1);
	skimTree->SetBranchAddress("pPAprimaryVertexFilter",&primaryVertexFilterBit,&primaryVertexBranch);
	skimTree->SetBranchStatus("pBeamScrapingFilter",1);
	skimTree->SetBranchAddress("pBeamScrapingFilter",&beamScrapingFilterBit,&beamScrapingBranch);
	skimTree->SetBranchStatus("phfCoincFilter",1);
	skimTree->SetBranchAddress("phfCoincFilter", &hfCoincidenceFilterBit, &hfCoincidenceBranch);
	skimTree->SetBranchStatus("pVertexFilterCutdz1p0",1);
	skimTree->SetBranchAddress("pVertexFilterCutdz1p0", &pVertexFilterCutdz1p0Bit, &pVertexFilterCutdz1p0Branch);
	skimTree->SetBranchStatus("pVertexFilterCutGplus",1);
	skimTree->SetBranchAddress("pVertexFilterCutGplus", &pVertexFilterCutGplusBit, &pVertexFilterCutGplusBranch);
	skimTree->SetBranchStatus("pVertexFilterCutVtx1",1);
	skimTree->SetBranchAddress("pVertexFilterCutVtx1", &pVertexFilterCutVtx1Bit, &pVertexFilterCutVtx1Branch);

	// Connect the branches to the track tree
	trackTree->SetBranchStatus("*",0);

	trackTree->SetBranchStatus("nTrk",1);
	trackTree->SetBranchAddress("nTrk",&nTracks,&nTracksBranch);
	trackTree->SetBranchStatus("highPurity",1);
	trackTree->SetBranchAddress("highPurity",&trackHighPurityArray,&trackHighPurityBranch);
	trackTree->SetBranchStatus("trkPt",1);
	trackTree->SetBranchAddress("trkPt",&trackPtArray,&trackPtBranch);
	trackTree->SetBranchStatus("trkPtError",1);
	trackTree->SetBranchAddress("trkPtError",&trackPtErrorArray,&trackPtErrorBranch);
	trackTree->SetBranchStatus("trkPhi",1);
	trackTree->SetBranchAddress("trkPhi",&trackPhiArray,&trackPhiBranch);
	trackTree->SetBranchStatus("trkEta",1);
	trackTree->SetBranchAddress("trkEta",&trackEtaArray,&trackEtaBranch);
	trackTree->SetBranchStatus("trkDz1",1);
	trackTree->SetBranchAddress("trkDz1",&trackVertexDistanceZArray,&trackVertexDistanceZBranch);
	trackTree->SetBranchStatus("trkDzError1",1);
	trackTree->SetBranchAddress("trkDzError1",&trackVertexDistanceZErrorArray,&trackVertexDistanceZErrorBranch);
	trackTree->SetBranchStatus("trkDxy1",1);
	trackTree->SetBranchAddress("trkDxy1",&trackVertexDistanceXYArray,&trackVertexDistanceXYBranch);
	trackTree->SetBranchStatus("trkDxyError1",1);
	trackTree->SetBranchAddress("trkDxyError1",&trackVertexDistanceXYErrorArray,&trackVertexDistanceXYErrorBranch);
	trackTree->SetBranchStatus("trkCharge",1);
	trackTree->SetBranchAddress("trkCharge",&trackChargeArray,&trackChargeBranch);


	// ========================================== //
	//			 Define output trees
	// ========================================== //


	//--------------- V0s ---------------
  	TTree *K0sTreeOutput = new TTree("K0sTree","");

	std::vector<double> *K0s_dxy1Vector = new std::vector<double>(); K0s_dxy1Vector->clear();
	std::vector<double> *K0s_dz1Vector= new std::vector<double>(); K0s_dz1Vector->clear();
	std::vector<double> *K0s_chi21Vector = new std::vector<double>(); K0s_chi21Vector->clear();
	std::vector<double> *K0s_d1pxVector= new std::vector<double>(); K0s_d1pxVector->clear();
	std::vector<double> *K0s_d1pyVector = new std::vector<double>(); K0s_d1pyVector->clear();
	std::vector<double> *K0s_d1pzVector= new std::vector<double>(); K0s_d1pzVector->clear();
	std::vector<double> *K0s_d1MVector = new std::vector<double>(); K0s_d1MVector->clear();
	std::vector<double> *K0s_d1NhitVector= new std::vector<double>(); K0s_d1NhitVector->clear();
	std::vector<double> *K0s_dxy2Vector = new std::vector<double>(); K0s_dxy2Vector->clear();
	std::vector<double> *K0s_dz2Vector= new std::vector<double>(); K0s_dz2Vector->clear();
	std::vector<double> *K0s_chi22Vector = new std::vector<double>(); K0s_chi22Vector->clear();
	std::vector<double> *K0s_d2pxVector= new std::vector<double>(); K0s_d2pxVector->clear();
	std::vector<double> *K0s_d2pyVector = new std::vector<double>(); K0s_d2pyVector->clear();
	std::vector<double> *K0s_d2pzVector= new std::vector<double>(); K0s_d2pzVector->clear();
	std::vector<double> *K0s_d2MVector = new std::vector<double>(); K0s_d2MVector->clear();
	std::vector<double> *K0s_d2NhitVector= new std::vector<double>(); K0s_d2NhitVector->clear();
	std::vector<double> *K0s_3DaglVector = new std::vector<double>(); K0s_3DaglVector->clear();
	std::vector<double> *K0s_3DdlVector= new std::vector<double>(); K0s_3DdlVector->clear();
	std::vector<double> *K0s_ptVector = new std::vector<double>(); K0s_ptVector->clear();
	std::vector<double> *K0s_etaVector= new std::vector<double>(); K0s_etaVector->clear();
	std::vector<double> *K0s_phiVector = new std::vector<double>(); K0s_phiVector->clear();
	std::vector<double> *K0s_massVector= new std::vector<double>(); K0s_massVector->clear();
	std::vector<double> *K0s_dcaVector = new std::vector<double>(); K0s_dcaVector->clear();
	std::vector<double> *K0s_vtxVector= new std::vector<double>(); K0s_vtxVector->clear();

	K0sTreeOutput->Branch("K0s_dxy1","vector<double>", &K0s_dxy1Vector);
	K0sTreeOutput->Branch("K0s_dz1","vector<double>", &K0s_dz1Vector);
	K0sTreeOutput->Branch("K0s_chi21","vector<double>", &K0s_chi21Vector);
	K0sTreeOutput->Branch("K0s_d1px","vector<double>", &K0s_d1pxVector);
	K0sTreeOutput->Branch("K0s_d1py","vector<double>", &K0s_d1pyVector);
	K0sTreeOutput->Branch("K0s_d1pz","vector<double>", &K0s_d1pzVector);
	K0sTreeOutput->Branch("K0s_d1M","vector<double>", &K0s_d1MVector);
	K0sTreeOutput->Branch("K0s_d1Nhit","vector<double>", &K0s_d1NhitVector);
	K0sTreeOutput->Branch("K0s_dxy2","vector<double>", &K0s_dxy2Vector);
	K0sTreeOutput->Branch("K0s_dz2","vector<double>", &K0s_dz2Vector);
	K0sTreeOutput->Branch("K0s_chi22","vector<double>", &K0s_chi22Vector);
	K0sTreeOutput->Branch("K0s_d2px","vector<double>", &K0s_d2pxVector);
	K0sTreeOutput->Branch("K0s_d2py","vector<double>", &K0s_d2pyVector);
	K0sTreeOutput->Branch("K0s_d2pz","vector<double>", &K0s_d2pzVector);
	K0sTreeOutput->Branch("K0s_d2M","vector<double>", &K0s_d2MVector);
	K0sTreeOutput->Branch("K0s_d2Nhit","vector<double>", &K0s_d2NhitVector);
	K0sTreeOutput->Branch("K0s_3Dagl","vector<double>", &K0s_3DaglVector);
	K0sTreeOutput->Branch("K0s_3Ddl","vector<double>", &K0s_3DdlVector);
	K0sTreeOutput->Branch("K0s_pt","vector<double>", &K0s_ptVector);
	K0sTreeOutput->Branch("K0s_eta","vector<double>", &K0s_etaVector);
	K0sTreeOutput->Branch("K0s_phi","vector<double>", &K0s_phiVector);
	K0sTreeOutput->Branch("K0s_mass","vector<double>", &K0s_massVector);
	K0sTreeOutput->Branch("K0s_dca","vector<double>", &K0s_dcaVector);
	K0sTreeOutput->Branch("K0s_vtx","vector<double>", &K0s_vtxVector);

  	TTree *LamTreeOutput = new TTree("LamTree","");

	std::vector<double> *Lam_dxy1Vector = new std::vector<double>(); Lam_dxy1Vector->clear();
	std::vector<double> *Lam_dz1Vector= new std::vector<double>(); Lam_dz1Vector->clear();
	std::vector<double> *Lam_chi21Vector = new std::vector<double>(); Lam_chi21Vector->clear();
	std::vector<double> *Lam_d1pxVector= new std::vector<double>(); Lam_d1pxVector->clear();
	std::vector<double> *Lam_d1pyVector = new std::vector<double>(); Lam_d1pyVector->clear();
	std::vector<double> *Lam_d1pzVector= new std::vector<double>(); Lam_d1pzVector->clear();
	std::vector<double> *Lam_d1MVector = new std::vector<double>(); Lam_d1MVector->clear();
	std::vector<double> *Lam_d1NhitVector= new std::vector<double>(); Lam_d1NhitVector->clear();

	std::vector<double> *Lam_dxy2Vector = new std::vector<double>(); Lam_dxy2Vector->clear();
	std::vector<double> *Lam_dz2Vector= new std::vector<double>(); Lam_dz2Vector->clear();
	std::vector<double> *Lam_chi22Vector = new std::vector<double>(); Lam_chi22Vector->clear();
	std::vector<double> *Lam_d2pxVector= new std::vector<double>(); Lam_d2pxVector->clear();
	std::vector<double> *Lam_d2pyVector = new std::vector<double>(); Lam_d2pyVector->clear();
	std::vector<double> *Lam_d2pzVector= new std::vector<double>(); Lam_d2pzVector->clear();
	std::vector<double> *Lam_d2MVector = new std::vector<double>(); Lam_d2MVector->clear();
	std::vector<double> *Lam_d2NhitVector= new std::vector<double>(); Lam_d2NhitVector->clear();

	std::vector<double> *Lam_3DaglVector = new std::vector<double>(); Lam_3DaglVector->clear();
	std::vector<double> *Lam_3DdlVector= new std::vector<double>(); Lam_3DdlVector->clear();
	std::vector<double> *Lam_ptVector = new std::vector<double>(); Lam_ptVector->clear();
	std::vector<double> *Lam_etaVector= new std::vector<double>(); Lam_etaVector->clear();
	std::vector<double> *Lam_phiVector = new std::vector<double>(); Lam_phiVector->clear();
	std::vector<double> *Lam_massVector= new std::vector<double>(); Lam_massVector->clear();
	std::vector<double> *Lam_dcaVector = new std::vector<double>(); Lam_dcaVector->clear();
	std::vector<double> *Lam_vtxVector= new std::vector<double>(); Lam_vtxVector->clear();
	std::vector<double> *Lam_idVector= new std::vector<double>(); Lam_idVector->clear();

	LamTreeOutput->Branch("Lam_dxy1","vector<double>", &Lam_dxy1Vector);
	LamTreeOutput->Branch("Lam_dz1","vector<double>", &Lam_dz1Vector);
	LamTreeOutput->Branch("Lam_chi21","vector<double>", &Lam_chi21Vector);
	LamTreeOutput->Branch("Lam_d1px","vector<double>", &Lam_d1pxVector);
	LamTreeOutput->Branch("Lam_d1py","vector<double>", &Lam_d1pyVector);
	LamTreeOutput->Branch("Lam_d1pz","vector<double>", &Lam_d1pzVector);
	LamTreeOutput->Branch("Lam_d1M","vector<double>", &Lam_d1MVector);
	LamTreeOutput->Branch("Lam_d1Nhit","vector<double>", &Lam_d1NhitVector);
	LamTreeOutput->Branch("Lam_dxy2","vector<double>", &Lam_dxy2Vector);
	LamTreeOutput->Branch("Lam_dz2","vector<double>", &Lam_dz2Vector);
	LamTreeOutput->Branch("Lam_chi22","vector<double>", &Lam_chi22Vector);
	LamTreeOutput->Branch("Lam_d2px","vector<double>", &Lam_d2pxVector);
	LamTreeOutput->Branch("Lam_d2py","vector<double>", &Lam_d2pyVector);
	LamTreeOutput->Branch("Lam_d2pz","vector<double>", &Lam_d2pzVector);
	LamTreeOutput->Branch("Lam_d2M","vector<double>", &Lam_d2MVector);
	LamTreeOutput->Branch("Lam_d2Nhit","vector<double>", &Lam_d2NhitVector);
	LamTreeOutput->Branch("Lam_3Dagl","vector<double>", &Lam_3DaglVector);
	LamTreeOutput->Branch("Lam_3Ddl","vector<double>", &Lam_3DdlVector);
	LamTreeOutput->Branch("Lam_pt","vector<double>", &Lam_ptVector);
	LamTreeOutput->Branch("Lam_eta","vector<double>", &Lam_etaVector);
	LamTreeOutput->Branch("Lam_phi","vector<double>", &Lam_phiVector);
	LamTreeOutput->Branch("Lam_mass","vector<double>", &Lam_massVector);
	LamTreeOutput->Branch("Lam_dca","vector<double>", &Lam_dcaVector);
	LamTreeOutput->Branch("Lam_vtx","vector<double>", &Lam_vtxVector);
	LamTreeOutput->Branch("Lam_id","vector<double>", &Lam_idVector);


   	TTree *CasTreeOutput = new TTree("XiTree","");

	std::vector<double> *Xi_d1ptVector = new std::vector<double>(); Xi_d1ptVector->clear();
	std::vector<double> *Xi_d1etaVector= new std::vector<double>(); Xi_d1etaVector->clear();
	std::vector<double> *Xi_d1phiVector = new std::vector<double>(); Xi_d1phiVector->clear();
	std::vector<double> *Xi_d1massVector = new std::vector<double>(); Xi_d1massVector->clear();

	std::vector<double> *Xi_chi21_1Vector = new std::vector<double>(); Xi_chi21_1Vector->clear();
	std::vector<double> *Xi_d1pt_1Vector= new std::vector<double>(); Xi_d1pt_1Vector->clear();
	std::vector<double> *Xi_d1eta_1Vector = new std::vector<double>(); Xi_d1eta_1Vector->clear();
	std::vector<double> *Xi_d1phi_1Vector = new std::vector<double>(); Xi_d1phi_1Vector->clear();
	std::vector<double> *Xi_d1mass_1Vector= new std::vector<double>(); Xi_d1mass_1Vector->clear();
	std::vector<double> *Xi_d1Nhit_1Vector = new std::vector<double>(); Xi_d1Nhit_1Vector->clear();

	std::vector<double> *Xi_chi21_2Vector= new std::vector<double>(); Xi_chi21_2Vector->clear();
	std::vector<double> *Xi_d1pt_2Vector = new std::vector<double>(); Xi_d1pt_2Vector->clear();
	std::vector<double> *Xi_d1eta_2Vector= new std::vector<double>(); Xi_d1eta_2Vector->clear();
	std::vector<double> *Xi_d1phi_2Vector = new std::vector<double>(); Xi_d1phi_2Vector->clear();
	std::vector<double> *Xi_d1mass_2Vector= new std::vector<double>(); Xi_d1mass_2Vector->clear();
	std::vector<double> *Xi_d1Nhit_2Vector = new std::vector<double>(); Xi_d1Nhit_2Vector->clear();

	std::vector<double> *Xi_chi22Vector= new std::vector<double>(); Xi_chi22Vector->clear();
	std::vector<double> *Xi_d2ptVector = new std::vector<double>(); Xi_d2ptVector->clear();
	std::vector<double> *Xi_d2etaVector= new std::vector<double>(); Xi_d2etaVector->clear();
	std::vector<double> *Xi_d2phiVector = new std::vector<double>(); Xi_d2phiVector->clear();
	std::vector<double> *Xi_d2massVector= new std::vector<double>(); Xi_d2massVector->clear();
	std::vector<double> *Xi_d2NhitVector = new std::vector<double>(); Xi_d2NhitVector->clear();

	std::vector<double> *Xi_cas3DIpSigValueVector = new std::vector<double>(); Xi_cas3DIpSigValueVector->clear();
	std::vector<double> *Xi_casPi3DIpSigValueVector = new std::vector<double>(); Xi_casPi3DIpSigValueVector->clear();
	std::vector<double> *Xi_VTrkPi3DIpSigValueVector = new std::vector<double>(); Xi_VTrkPi3DIpSigValueVector->clear();
	std::vector<double> *Xi_VTrkP3DIpSigValueVector = new std::vector<double>(); Xi_VTrkP3DIpSigValueVector->clear();
	std::vector<double> *Xi_casFlightSigValueVector = new std::vector<double>(); Xi_casFlightSigValueVector->clear();
	std::vector<double> *Xi_distanceSigValueVector= new std::vector<double>(); Xi_distanceSigValueVector->clear();
	std::vector<double> *Xi_ptVector = new std::vector<double>(); Xi_ptVector->clear();
	std::vector<double> *Xi_etaVector= new std::vector<double>(); Xi_etaVector->clear();
	std::vector<double> *Xi_phiVector = new std::vector<double>(); Xi_phiVector->clear();
	std::vector<double> *Xi_massVector= new std::vector<double>(); Xi_massVector->clear();
	std::vector<double> *Xi_idVector= new std::vector<double>(); Xi_idVector->clear();

	CasTreeOutput->Branch("Xi_d1pt","vector<double>", &Xi_d1ptVector);
	CasTreeOutput->Branch("Xi_d1eta","vector<double>", &Xi_d1etaVector);
	CasTreeOutput->Branch("Xi_d1phi","vector<double>", &Xi_d1phiVector);
	CasTreeOutput->Branch("Xi_d1mass","vector<double>", &Xi_d1massVector);
	CasTreeOutput->Branch("Xi_chi21_1","vector<double>", &Xi_chi21_1Vector);
	CasTreeOutput->Branch("Xi_d1pt_1","vector<double>", &Xi_d1pt_1Vector);
	CasTreeOutput->Branch("Xi_d1eta_1","vector<double>", &Xi_d1eta_1Vector);
	CasTreeOutput->Branch("Xi_d1phi_1","vector<double>", &Xi_d1phi_1Vector);
	CasTreeOutput->Branch("Xi_d1mass_1","vector<double>", &Xi_d1mass_1Vector);
	CasTreeOutput->Branch("Xi_d1Nhit_1","vector<double>", &Xi_d1Nhit_1Vector);
	CasTreeOutput->Branch("Xi_chi21_2","vector<double>", &Xi_chi21_2Vector);
	CasTreeOutput->Branch("Xi_d1pt_2","vector<double>", &Xi_d1pt_2Vector);
	CasTreeOutput->Branch("Xi_d1eta_2","vector<double>", &Xi_d1eta_2Vector);
	CasTreeOutput->Branch("Xi_d1phi_2","vector<double>", &Xi_d1phi_2Vector);
	CasTreeOutput->Branch("Xi_d1mass_2","vector<double>", &Xi_d1mass_2Vector);
	CasTreeOutput->Branch("Xi_d1Nhit_2","vector<double>", &Xi_d1Nhit_2Vector);
	CasTreeOutput->Branch("Xi_chi22","vector<double>", &Xi_chi22Vector);
	CasTreeOutput->Branch("Xi_d2pt","vector<double>", &Xi_d2ptVector);
	CasTreeOutput->Branch("Xi_d2eta","vector<double>", &Xi_d2etaVector);
	CasTreeOutput->Branch("Xi_d2phi","vector<double>", &Xi_d2phiVector);
	CasTreeOutput->Branch("Xi_d2mass","vector<double>", &Xi_d2massVector);
	CasTreeOutput->Branch("Xi_d2Nhit","vector<double>", &Xi_d2NhitVector);
	CasTreeOutput->Branch("Xi_cas3DIpSigValue","vector<double>", &Xi_cas3DIpSigValueVector);
	CasTreeOutput->Branch("Xi_casPi3DIpSigValue","vector<double>", &Xi_casPi3DIpSigValueVector);
	CasTreeOutput->Branch("Xi_VTrkPi3DIpSigValue","vector<double>", &Xi_VTrkPi3DIpSigValueVector);
	CasTreeOutput->Branch("Xi_VTrkP3DIpSigValue","vector<double>", &Xi_VTrkP3DIpSigValueVector);
	CasTreeOutput->Branch("Xi_casFlightSigValue","vector<double>", &Xi_casFlightSigValueVector);
	CasTreeOutput->Branch("Xi_distanceSigValue","vector<double>", &Xi_distanceSigValueVector);
	CasTreeOutput->Branch("Xi_pt","vector<double>", &Xi_ptVector);
	CasTreeOutput->Branch("Xi_eta","vector<double>", &Xi_etaVector);
	CasTreeOutput->Branch("Xi_phi","vector<double>", &Xi_phiVector);
	CasTreeOutput->Branch("Xi_mass","vector<double>", &Xi_massVector);
	CasTreeOutput->Branch("Xi_id","vector<double>", &Xi_idVector);


   	TTree *OmeTreeOutput = new TTree("OmTree","");

	std::vector<double> *Om_d1ptVector = new std::vector<double>(); Om_d1ptVector->clear();
	std::vector<double> *Om_d1etaVector= new std::vector<double>(); Om_d1etaVector->clear();
	std::vector<double> *Om_d1phiVector = new std::vector<double>(); Om_d1phiVector->clear();
	std::vector<double> *Om_d1massVector = new std::vector<double>(); Om_d1massVector->clear();

	std::vector<double> *Om_chi21_1Vector = new std::vector<double>(); Om_chi21_1Vector->clear();
	std::vector<double> *Om_d1pt_1Vector= new std::vector<double>(); Om_d1pt_1Vector->clear();
	std::vector<double> *Om_d1eta_1Vector = new std::vector<double>(); Om_d1eta_1Vector->clear();
	std::vector<double> *Om_d1phi_1Vector = new std::vector<double>(); Om_d1phi_1Vector->clear();
	std::vector<double> *Om_d1mass_1Vector= new std::vector<double>(); Om_d1mass_1Vector->clear();
	std::vector<double> *Om_d1Nhit_1Vector = new std::vector<double>(); Om_d1Nhit_1Vector->clear();

	std::vector<double> *Om_chi21_2Vector= new std::vector<double>(); Om_chi21_2Vector->clear();
	std::vector<double> *Om_d1pt_2Vector = new std::vector<double>(); Om_d1pt_2Vector->clear();
	std::vector<double> *Om_d1eta_2Vector= new std::vector<double>(); Om_d1eta_2Vector->clear();
	std::vector<double> *Om_d1phi_2Vector = new std::vector<double>(); Om_d1phi_2Vector->clear();
	std::vector<double> *Om_d1mass_2Vector= new std::vector<double>(); Om_d1mass_2Vector->clear();
	std::vector<double> *Om_d1Nhit_2Vector = new std::vector<double>(); Om_d1Nhit_2Vector->clear();

	std::vector<double> *Om_chi22Vector= new std::vector<double>(); Om_chi22Vector->clear();
	std::vector<double> *Om_d2ptVector = new std::vector<double>(); Om_d2ptVector->clear();
	std::vector<double> *Om_d2etaVector= new std::vector<double>(); Om_d2etaVector->clear();
	std::vector<double> *Om_d2phiVector = new std::vector<double>(); Om_d2phiVector->clear();
	std::vector<double> *Om_d2massVector= new std::vector<double>(); Om_d2massVector->clear();
	std::vector<double> *Om_d2NhitVector = new std::vector<double>(); Om_d2NhitVector->clear();

	std::vector<double> *Om_cas3DIpSigValueVector = new std::vector<double>(); Om_cas3DIpSigValueVector->clear();
	std::vector<double> *Om_casPi3DIpSigValueVector = new std::vector<double>(); Om_casPi3DIpSigValueVector->clear();
	std::vector<double> *Om_VTrkPi3DIpSigValueVector = new std::vector<double>(); Om_VTrkPi3DIpSigValueVector->clear();
	std::vector<double> *Om_VTrkP3DIpSigValueVector = new std::vector<double>(); Om_VTrkP3DIpSigValueVector->clear();
	std::vector<double> *Om_casFlightSigValueVector = new std::vector<double>(); Om_casFlightSigValueVector->clear();
	std::vector<double> *Om_distanceSigValueVector= new std::vector<double>(); Om_distanceSigValueVector->clear();
	std::vector<double> *Om_ptVector = new std::vector<double>(); Om_ptVector->clear();
	std::vector<double> *Om_etaVector= new std::vector<double>(); Om_etaVector->clear();
	std::vector<double> *Om_phiVector = new std::vector<double>(); Om_phiVector->clear();
	std::vector<double> *Om_massVector= new std::vector<double>(); Om_massVector->clear();
	std::vector<double> *Om_idVector= new std::vector<double>(); Om_idVector->clear();

	OmeTreeOutput->Branch("Om_d1pt","vector<double>", &Om_d1ptVector);
	OmeTreeOutput->Branch("Om_d1eta","vector<double>", &Om_d1etaVector);
	OmeTreeOutput->Branch("Om_d1phi","vector<double>", &Om_d1phiVector);
	OmeTreeOutput->Branch("Om_d1mass","vector<double>", &Om_d1massVector);
	OmeTreeOutput->Branch("Om_chi21_1","vector<double>", &Om_chi21_1Vector);
	OmeTreeOutput->Branch("Om_d1pt_1","vector<double>", &Om_d1pt_1Vector);
	OmeTreeOutput->Branch("Om_d1eta_1","vector<double>", &Om_d1eta_1Vector);
	OmeTreeOutput->Branch("Om_d1phi_1","vector<double>", &Om_d1phi_1Vector);
	OmeTreeOutput->Branch("Om_d1mass_1","vector<double>", &Om_d1mass_1Vector);
	OmeTreeOutput->Branch("Om_d1Nhit_1","vector<double>", &Om_d1Nhit_1Vector);
	OmeTreeOutput->Branch("Om_chi21_2","vector<double>", &Om_chi21_2Vector);
	OmeTreeOutput->Branch("Om_d1pt_2","vector<double>", &Om_d1pt_2Vector);
	OmeTreeOutput->Branch("Om_d1eta_2","vector<double>", &Om_d1eta_2Vector);
	OmeTreeOutput->Branch("Om_d1phi_2","vector<double>", &Om_d1phi_2Vector);
	OmeTreeOutput->Branch("Om_d1mass_2","vector<double>", &Om_d1mass_2Vector);
	OmeTreeOutput->Branch("Om_d1Nhit_2","vector<double>", &Om_d1Nhit_2Vector);
	OmeTreeOutput->Branch("Om_chi22","vector<double>", &Om_chi22Vector);
	OmeTreeOutput->Branch("Om_d2pt","vector<double>", &Om_d2ptVector);
	OmeTreeOutput->Branch("Om_d2eta","vector<double>", &Om_d2etaVector);
	OmeTreeOutput->Branch("Om_d2phi","vector<double>", &Om_d2phiVector);
	OmeTreeOutput->Branch("Om_d2mass","vector<double>", &Om_d2massVector);
	OmeTreeOutput->Branch("Om_d2Nhit","vector<double>", &Om_d2NhitVector);
	OmeTreeOutput->Branch("Om_cas3DIpSigValue","vector<double>", &Om_cas3DIpSigValueVector);
	OmeTreeOutput->Branch("Om_casPi3DIpSigValue","vector<double>", &Om_casPi3DIpSigValueVector);
	OmeTreeOutput->Branch("Om_VTrkPi3DIpSigValue","vector<double>", &Om_VTrkPi3DIpSigValueVector);
	OmeTreeOutput->Branch("Om_VTrkP3DIpSigValue","vector<double>", &Om_VTrkP3DIpSigValueVector);
	OmeTreeOutput->Branch("Om_casFlightSigValue","vector<double>", &Om_casFlightSigValueVector);
	OmeTreeOutput->Branch("Om_distanceSigValue","vector<double>", &Om_distanceSigValueVector);
	OmeTreeOutput->Branch("Om_pt","vector<double>", &Om_ptVector);
	OmeTreeOutput->Branch("Om_eta","vector<double>", &Om_etaVector);
	OmeTreeOutput->Branch("Om_phi","vector<double>", &Om_phiVector);
	OmeTreeOutput->Branch("Om_mass","vector<double>", &Om_massVector);
	OmeTreeOutput->Branch("Om_id","vector<double>", &Om_idVector);


  
	// Copy the heavy ion tree to the output
	TTree *heavyIonTreeOutput = new TTree("HiTree","");
	// Connect the branches of the heavy ion tree
	heavyIonTreeOutput->Branch("run",&run,"run/i");
	heavyIonTreeOutput->Branch("evt",&event,"evt/l");
	heavyIonTreeOutput->Branch("lumi",&lumi,"lumi/i");
	heavyIonTreeOutput->Branch("vz",&vertexZ,"vz/F");
	heavyIonTreeOutput->Branch("hiHFplus",&hiHFplus,"hiHFplus/F");
	heavyIonTreeOutput->Branch("hiHFminus",&hiHFminus,"hiHFminus/F");
	heavyIonTreeOutput->Branch("hiZDCplus",&hiZDCplus,"hiZDCplus/F");
	heavyIonTreeOutput->Branch("hiZDCminus",&hiZDCminus,"hiZDCminus/F");
	heavyIonTreeOutput->Branch("hi_FRG",&hi_FRG,"hi_FRG/F");
	heavyIonTreeOutput->Branch("Ntrkoff",&Ntrkoff,"Ntrkoff/I");
	heavyIonTreeOutput->Branch("NtrkCorr",&NtrkCorr,"NtrkCorr/F");
	heavyIonTreeOutput->Branch("NtrkCorrTight",&NtrkCorrTight,"NtrkCorrTight/F");
	heavyIonTreeOutput->Branch("NtrkCorrLoose",&NtrkCorrLoose,"NtrkCorrLoose/F");

	// Copy the skim tree to the output
	TTree *skimTreeOutput = new TTree("HltTree","");
	skimTreeOutput->Branch("pPAprimaryVertexFilter",&primaryVertexFilterBit,"pPAprimaryVertexFilter/I");
	skimTreeOutput->Branch("pBeamScrapingFilter",&beamScrapingFilterBit,"pBeamScrapingFilter/I");
	skimTreeOutput->Branch("phfCoincFilter", &hfCoincidenceFilterBit, "phfCoincFilter/I");
	skimTreeOutput->Branch("pVertexFilterCutdz1p0", &pVertexFilterCutdz1p0Bit, "pVertexFilterCutdz1p0/I");
	skimTreeOutput->Branch("pVertexFilterCutGplus",&pVertexFilterCutGplusBit,"pVertexFilterCutGplus/I");
	skimTreeOutput->Branch("pVertexFilterCutVtx1",&pVertexFilterCutVtx1Bit,"pVertexFilterCutVtx1/I");

	// reading the track efficiency files
	TFile *fileeff_nominal = TFile::Open("Hijing_8TeV_dataBS.root","READ");
    TH2 *eff_factor_nominal = nullptr; 
    fileeff_nominal->GetObject("rTotalEff3D_0", eff_factor_nominal);  // eff / (1 - fake rate)

	TFile *fileeff_tight = TFile::Open("Hijing_8TeV_MB_eff_v3_tight.root","READ");
    TH2 *eff_factor_tight = nullptr; 
    fileeff_tight->GetObject("rTotalEff3D_0", eff_factor_tight);  // eff / (1 - fake rate)

	TFile *fileeff_loose = TFile::Open("Hijing_8TeV_MB_eff_v3_loose.root","READ");
    TH2 *eff_factor_loose = nullptr; 
    fileeff_loose->GetObject("rTotalEff3D_0", eff_factor_loose);  // eff / (1 - fake rate)



	// ========================================== //
	//			Starting matching events	      //
	// ========================================== //

	Int_t jettrk_events = heavyIonTree->GetEntries(); // number of events
    // loop through jets and create a key for each event
    for(int i_entry = 0; i_entry < jettrk_events; i_entry++){
       heavyIonTree->GetEntry(i_entry);
       unsigned long long key = keyFromRunLumiEvent(run, lumi, event);
       runLumiEvtToEntryMap[key] = i_entry;
    }

	// ========================================== //
	//				Loop over all events 		  //
	// ========================================== //

	int nEvents = MainV0Tree->GetEntries();
	cout << "There are " << nEvents << " events" << endl;

	for(int iEvent = 0; iEvent < nEvents; iEvent++) {
		
		if( iEvent % 1000 == 0 )	std::cout << "iEvent: " << iEvent <<	" of " << nEvents << std::endl;

		// ========================================== //
		//			Start with the V0s	              //
		// ========================================== //
		MainV0Tree->GetEntry(iEvent);
		K0sTree->GetEntry(iEvent);
		LamTree->GetEntry(iEvent);
		CasTree->GetEntry(iEvent);
		OmeTree->GetEntry(iEvent);

		//Find matching jet event
		if (V0_evt < 0) continue;
		if (V0_evt > INT_MAX - 1) continue;
		unsigned long long key = keyFromRunLumiEvent((UInt_t)V0_run,(UInt_t)V0_lumi,(ULong64_t)V0_evt);
       	long long i_entry = -1;
       	if(runLumiEvtToEntryMap.count(key) == 0) continue; // skip reco event if there is no event match
       	else i_entry = runLumiEvtToEntryMap.at(key);
		
		// ========================================== //
		//	Read the event to input trees	      //
		// ========================================== //
		
		heavyIonTree->GetEntry(i_entry);
		skimTree->GetEntry(i_entry);
		trackTree->GetEntry(i_entry);

		Ntrkoff = get_Ntrkoff(nTracks, trackEtaArray, trackPtArray, trackChargeArray, trackHighPurityArray, trackPtErrorArray, trackVertexDistanceXYArray, trackVertexDistanceXYErrorArray, trackVertexDistanceZArray, trackVertexDistanceZErrorArray);
		NtrkCorr = get_Ntrkcorr(eff_factor_nominal, nTracks, trackEtaArray, trackPtArray, trackChargeArray, trackHighPurityArray, trackPtErrorArray, trackVertexDistanceXYArray, trackVertexDistanceXYErrorArray, trackVertexDistanceZArray, trackVertexDistanceZErrorArray);
		NtrkCorrTight = get_Ntrkcorr_tight(eff_factor_tight, nTracks, trackEtaArray, trackPtArray, trackChargeArray, trackHighPurityArray, trackPtErrorArray, trackVertexDistanceXYArray, trackVertexDistanceXYErrorArray, trackVertexDistanceZArray, trackVertexDistanceZErrorArray);
		NtrkCorrLoose = get_Ntrkcorr_loose(eff_factor_loose, nTracks, trackEtaArray, trackPtArray, trackChargeArray, trackHighPurityArray, trackPtErrorArray, trackVertexDistanceXYArray, trackVertexDistanceXYErrorArray, trackVertexDistanceZArray, trackVertexDistanceZErrorArray);

		heavyIonTreeOutput->Fill(); // fill event information
		hltTreeOutput->Fill();      // HLT information
		skimTreeOutput->Fill();		// filter information

   		//loop over K0s
   		for(int iK0s = 0; iK0s < K0s_pt->size(); iK0s++){

			if(fabs(K0s_eta->at(iK0s)) > V0etamin) continue; //eta acceptance
			if(K0s_pt->at(iK0s) < V0ptmin) continue;   //Minimum V0 pT

			K0s_dxy1Vector->push_back(K0s_dxy1->at(iK0s));
			K0s_dz1Vector->push_back(K0s_dz1->at(iK0s));
			K0s_chi21Vector->push_back(K0s_chi21->at(iK0s));
			K0s_d1pxVector->push_back(K0s_d1px->at(iK0s));
			K0s_d1pyVector->push_back(K0s_d1py->at(iK0s));
			K0s_d1pzVector->push_back(K0s_d1pz->at(iK0s));
			K0s_d1MVector->push_back(K0s_d1M->at(iK0s));
			K0s_d1NhitVector->push_back(K0s_d1Nhit->at(iK0s));

			K0s_dxy2Vector->push_back(K0s_dxy2->at(iK0s));
			K0s_dz2Vector->push_back(K0s_dz2->at(iK0s));
			K0s_chi22Vector->push_back(K0s_chi22->at(iK0s));
			K0s_d2pxVector->push_back(K0s_d2px->at(iK0s));
			K0s_d2pyVector->push_back(K0s_d2py->at(iK0s));
			K0s_d2pzVector->push_back(K0s_d2pz->at(iK0s));
			K0s_d2MVector->push_back(K0s_d2M->at(iK0s));
			K0s_d2NhitVector->push_back(K0s_d2Nhit->at(iK0s));

			K0s_3DaglVector->push_back(K0s_3Dagl->at(iK0s));
			K0s_3DdlVector->push_back(K0s_3Ddl->at(iK0s));
			K0s_ptVector->push_back(K0s_pt->at(iK0s));
			K0s_etaVector->push_back(K0s_eta->at(iK0s));
			K0s_phiVector->push_back(K0s_phi->at(iK0s));
			K0s_massVector->push_back(K0s_mass->at(iK0s));
			K0s_dcaVector->push_back(K0s_dca->at(iK0s));
			K0s_vtxVector->push_back(K0s_vtx->at(iK0s));

      		}
      
      	K0sTreeOutput->Fill();

    		// Clear the vectors before the next event! Otherwise all the K0s pile up cumulatively
		K0s_dxy1Vector->clear();
		K0s_dz1Vector->clear();
		K0s_chi21Vector->clear();
		K0s_d1pxVector->clear();
		K0s_d1pyVector->clear();
		K0s_d1pzVector->clear();
		K0s_d1MVector->clear();
		K0s_d1NhitVector->clear();

		K0s_dxy2Vector->clear();
		K0s_dz2Vector->clear();
		K0s_chi22Vector->clear();
		K0s_d2pxVector->clear();
		K0s_d2pyVector->clear();
		K0s_d2pzVector->clear();
		K0s_d2MVector->clear();
		K0s_d2NhitVector->clear();

		K0s_3DaglVector->clear();
		K0s_3DdlVector->clear();
		K0s_ptVector->clear();
		K0s_etaVector->clear();
		K0s_phiVector->clear();
		K0s_massVector->clear();
		K0s_dcaVector->clear();
		K0s_vtxVector->clear();
 
    	//loop over Lambdas
   		for(int iLam = 0; iLam < Lam_pt->size(); iLam++){

			if(fabs(Lam_eta->at(iLam)) > V0etamin) continue; //eta acceptance
			if(Lam_pt->at(iLam) < V0ptmin) continue;   //Minimum V0 pT

			Lam_dxy1Vector->push_back(Lam_dxy1->at(iLam));
			Lam_dz1Vector->push_back(Lam_dz1->at(iLam));
			Lam_chi21Vector->push_back(Lam_chi21->at(iLam));
			Lam_d1pxVector->push_back(Lam_d1px->at(iLam));
			Lam_d1pyVector->push_back(Lam_d1py->at(iLam));
			Lam_d1pzVector->push_back(Lam_d1pz->at(iLam));
			Lam_d1MVector->push_back(Lam_d1M->at(iLam));
			Lam_d1NhitVector->push_back(Lam_d1Nhit->at(iLam));

			Lam_dxy2Vector->push_back(Lam_dxy2->at(iLam));
			Lam_dz2Vector->push_back(Lam_dz2->at(iLam));
			Lam_chi22Vector->push_back(Lam_chi22->at(iLam));
			Lam_d2pxVector->push_back(Lam_d2px->at(iLam));
			Lam_d2pyVector->push_back(Lam_d2py->at(iLam));
			Lam_d2pzVector->push_back(Lam_d2pz->at(iLam));
			Lam_d2MVector->push_back(Lam_d2M->at(iLam));
			Lam_d2NhitVector->push_back(Lam_d2Nhit->at(iLam));

			Lam_3DaglVector->push_back(Lam_3Dagl->at(iLam));
			Lam_3DdlVector->push_back(Lam_3Ddl->at(iLam));
			Lam_ptVector->push_back(Lam_pt->at(iLam));
			Lam_etaVector->push_back(Lam_eta->at(iLam));
			Lam_phiVector->push_back(Lam_phi->at(iLam));
			Lam_massVector->push_back(Lam_mass->at(iLam));
			Lam_dcaVector->push_back(Lam_dca->at(iLam));
			Lam_vtxVector->push_back(Lam_vtx->at(iLam));
			Lam_idVector->push_back(Lam_id->at(iLam));

      	}
      
   		LamTreeOutput->Fill();

    		// Clear the vectors before the next event! Otherwise all the Lam pile up cumulatively
		Lam_dxy1Vector->clear();
		Lam_dz1Vector->clear();
		Lam_chi21Vector->clear();
		Lam_d1pxVector->clear();
		Lam_d1pyVector->clear();
		Lam_d1pzVector->clear();
		Lam_d1MVector->clear();
		Lam_d1NhitVector->clear();

		Lam_dxy2Vector->clear();
		Lam_dz2Vector->clear();
		Lam_chi22Vector->clear();
		Lam_d2pxVector->clear();
		Lam_d2pyVector->clear();
		Lam_d2pzVector->clear();
		Lam_d2MVector->clear();
		Lam_d2NhitVector->clear();

		Lam_3DaglVector->clear();
		Lam_3DdlVector->clear();
		Lam_ptVector->clear();
		Lam_etaVector->clear();
		Lam_phiVector->clear();
		Lam_massVector->clear();
		Lam_dcaVector->clear();
		Lam_vtxVector->clear();
		Lam_idVector->clear();

    	//loop over Xi
   		for(int iXi = 0; iXi < Xi_pt->size(); iXi++){

			if(fabs(Xi_eta->at(iXi)) > V0etamin) continue; //eta acceptance
			if(Xi_pt->at(iXi) < V0ptmin) continue;   //Minimum V0 pT

			Xi_d1ptVector->push_back(Xi_d1pt->at(iXi));
			Xi_d1etaVector->push_back(Xi_d1eta->at(iXi));
			Xi_d1phiVector->push_back(Xi_d1phi->at(iXi));
			Xi_d1massVector->push_back(Xi_d1mass->at(iXi));

			Xi_chi21_1Vector->push_back(Xi_chi21_1->at(iXi));
			Xi_d1pt_1Vector->push_back(Xi_d1pt_1->at(iXi));
			Xi_d1eta_1Vector->push_back(Xi_d1eta_1->at(iXi));
			Xi_d1phi_1Vector->push_back(Xi_d1phi_1->at(iXi));
			Xi_d1mass_1Vector->push_back(Xi_d1mass_1->at(iXi));
			Xi_d1Nhit_1Vector->push_back(Xi_d1Nhit_1->at(iXi));

			Xi_chi21_2Vector->push_back(Xi_chi21_2->at(iXi));
			Xi_d1pt_2Vector->push_back(Xi_d1pt_2->at(iXi));
			Xi_d1eta_2Vector->push_back(Xi_d1eta_2->at(iXi));
			Xi_d1phi_2Vector->push_back(Xi_d1phi_2->at(iXi));
			Xi_d1mass_2Vector->push_back(Xi_d1mass_2->at(iXi));
			Xi_d1Nhit_2Vector->push_back(Xi_d1Nhit_2->at(iXi));

			Xi_chi22Vector->push_back(Xi_chi22->at(iXi));
			Xi_d2ptVector->push_back(Xi_d2pt->at(iXi));
			Xi_d2etaVector->push_back(Xi_d2eta->at(iXi));
			Xi_d2phiVector->push_back(Xi_d2phi->at(iXi));
			Xi_d2massVector->push_back(Xi_d2mass->at(iXi));
			Xi_d2NhitVector->push_back(Xi_d2Nhit->at(iXi));

			Xi_cas3DIpSigValueVector->push_back(Xi_cas3DIpSigValue->at(iXi));
			Xi_casPi3DIpSigValueVector->push_back(Xi_casPi3DIpSigValue->at(iXi));
			Xi_VTrkPi3DIpSigValueVector->push_back(Xi_VTrkPi3DIpSigValue->at(iXi));
			Xi_VTrkP3DIpSigValueVector->push_back(Xi_VTrkP3DIpSigValue->at(iXi));
			Xi_casFlightSigValueVector->push_back(Xi_casFlightSigValue->at(iXi));
			Xi_distanceSigValueVector->push_back(Xi_distanceSigValue->at(iXi));
			Xi_ptVector->push_back(Xi_pt->at(iXi));
			Xi_etaVector->push_back(Xi_eta->at(iXi));
			Xi_phiVector->push_back(Xi_phi->at(iXi));
			Xi_massVector->push_back(Xi_mass->at(iXi));
			Xi_idVector->push_back(Xi_id->at(iXi));

      	}
      
   		CasTreeOutput->Fill();

    		// Clear the vectors before the next event! Otherwise all the Xi pile up cumulatively
		Xi_d1ptVector->clear();
		Xi_d1etaVector->clear();
		Xi_d1phiVector->clear();
		Xi_d1massVector->clear();

		Xi_chi21_1Vector->clear();
		Xi_d1pt_1Vector->clear();
		Xi_d1eta_1Vector->clear();
		Xi_d1phi_1Vector->clear();
		Xi_d1mass_1Vector->clear();
		Xi_d1Nhit_1Vector->clear();

		Xi_chi21_2Vector->clear();
		Xi_d1pt_2Vector->clear();
		Xi_d1eta_2Vector->clear();
		Xi_d1phi_2Vector->clear();
		Xi_d1mass_2Vector->clear();
		Xi_d1Nhit_2Vector->clear();

		Xi_chi22Vector->clear();
		Xi_d2ptVector->clear();
		Xi_d2etaVector->clear();
		Xi_d2phiVector->clear();
		Xi_d2massVector->clear();
		Xi_d2NhitVector->clear();

		Xi_cas3DIpSigValueVector->clear();
		Xi_casPi3DIpSigValueVector->clear();
		Xi_VTrkPi3DIpSigValueVector->clear();
		Xi_VTrkP3DIpSigValueVector->clear();
		Xi_casFlightSigValueVector->clear();
		Xi_distanceSigValueVector->clear();
		Xi_ptVector->clear();
		Xi_etaVector->clear();
		Xi_phiVector->clear();
		Xi_massVector->clear();
		Xi_idVector->clear();

    	//loop over Omegas
   		for(int iOm = 0; iOm < Om_pt->size(); iOm++){

			if(fabs(Om_eta->at(iOm)) > V0etamin) continue; //eta acceptance
			if(Om_pt->at(iOm) < V0ptmin) continue;   //Minimum V0 pT

			Om_d1ptVector->push_back(Om_d1pt->at(iOm));
			Om_d1etaVector->push_back(Om_d1eta->at(iOm));
			Om_d1phiVector->push_back(Om_d1phi->at(iOm));
			Om_d1massVector->push_back(Om_d1mass->at(iOm));

			Om_chi21_1Vector->push_back(Om_chi21_1->at(iOm));
			Om_d1pt_1Vector->push_back(Om_d1pt_1->at(iOm));
			Om_d1eta_1Vector->push_back(Om_d1eta_1->at(iOm));
			Om_d1phi_1Vector->push_back(Om_d1phi_1->at(iOm));
			Om_d1mass_1Vector->push_back(Om_d1mass_1->at(iOm));
			Om_d1Nhit_1Vector->push_back(Om_d1Nhit_1->at(iOm));

			Om_chi21_2Vector->push_back(Om_chi21_2->at(iOm));
			Om_d1pt_2Vector->push_back(Om_d1pt_2->at(iOm));
			Om_d1eta_2Vector->push_back(Om_d1eta_2->at(iOm));
			Om_d1phi_2Vector->push_back(Om_d1phi_2->at(iOm));
			Om_d1mass_2Vector->push_back(Om_d1mass_2->at(iOm));
			Om_d1Nhit_2Vector->push_back(Om_d1Nhit_2->at(iOm));

			Om_chi22Vector->push_back(Om_chi22->at(iOm));
			Om_d2ptVector->push_back(Om_d2pt->at(iOm));
			Om_d2etaVector->push_back(Om_d2eta->at(iOm));
			Om_d2phiVector->push_back(Om_d2phi->at(iOm));
			Om_d2massVector->push_back(Om_d2mass->at(iOm));
			Om_d2NhitVector->push_back(Om_d2Nhit->at(iOm));

			Om_cas3DIpSigValueVector->push_back(Om_cas3DIpSigValue->at(iOm));
			Om_casPi3DIpSigValueVector->push_back(Om_casPi3DIpSigValue->at(iOm));
			Om_VTrkPi3DIpSigValueVector->push_back(Om_VTrkPi3DIpSigValue->at(iOm));
			Om_VTrkP3DIpSigValueVector->push_back(Om_VTrkP3DIpSigValue->at(iOm));
			Om_casFlightSigValueVector->push_back(Om_casFlightSigValue->at(iOm));
			Om_distanceSigValueVector->push_back(Om_distanceSigValue->at(iOm));
			Om_ptVector->push_back(Om_pt->at(iOm));
			Om_etaVector->push_back(Om_eta->at(iOm));
			Om_phiVector->push_back(Om_phi->at(iOm));
			Om_massVector->push_back(Om_mass->at(iOm));
			Om_idVector->push_back(Om_id->at(iOm));

      	}
      
   		OmeTreeOutput->Fill();

    		// Clear the vectors before the next event! Otherwise all the Om pile up cumulatively
		Om_d1ptVector->clear();
		Om_d1etaVector->clear();
		Om_d1phiVector->clear();
		Om_d1massVector->clear();

		Om_chi21_1Vector->clear();
		Om_d1pt_1Vector->clear();
		Om_d1eta_1Vector->clear();
		Om_d1phi_1Vector->clear();
		Om_d1mass_1Vector->clear();
		Om_d1Nhit_1Vector->clear();

		Om_chi21_2Vector->clear();
		Om_d1pt_2Vector->clear();
		Om_d1eta_2Vector->clear();
		Om_d1phi_2Vector->clear();
		Om_d1mass_2Vector->clear();
		Om_d1Nhit_2Vector->clear();

		Om_chi22Vector->clear();
		Om_d2ptVector->clear();
		Om_d2etaVector->clear();
		Om_d2phiVector->clear();
		Om_d2massVector->clear();
		Om_d2NhitVector->clear();

		Om_cas3DIpSigValueVector->clear();
		Om_casPi3DIpSigValueVector->clear();
		Om_VTrkPi3DIpSigValueVector->clear();
		Om_VTrkP3DIpSigValueVector->clear();
		Om_casFlightSigValueVector->clear();
		Om_distanceSigValueVector->clear();
		Om_ptVector->clear();
		Om_etaVector->clear();
		Om_phiVector->clear();
		Om_massVector->clear();
		Om_idVector->clear();

	} // End loop over events

	// Write the skimmed trees to the output file
  	TFile *outputFile = new TFile(outputFileName, "RECREATE");

	gDirectory->mkdir("hiEvtAnalyzer");
	gDirectory->cd("hiEvtAnalyzer");
	heavyIonTreeOutput->Write();
	gDirectory->cd("../");

	gDirectory->mkdir("skimanalysis");
	gDirectory->cd("skimanalysis");
	skimTreeOutput->Write();
	gDirectory->cd("../");

	gDirectory->mkdir("K0sTree");
	gDirectory->cd("K0sTree");
	K0sTreeOutput->Write();	
	gDirectory->cd("../");

	gDirectory->mkdir("LamTree");
	gDirectory->cd("LamTree");
	LamTreeOutput->Write();	
	gDirectory->cd("../");

 	gDirectory->mkdir("XiTree");
	gDirectory->cd("XiTree");
	CasTreeOutput->Write();	
	gDirectory->cd("../");

 	gDirectory->mkdir("OmTree");
	gDirectory->cd("OmTree");
	OmeTreeOutput->Write();	
	gDirectory->cd("../");

	outputFile->Close();

	cout << endl;
	cout << "------------------------------------- SKIMMER DONE --------------------------------------" << endl;
	cout << endl;


	sec_end = clock(); // stop time counting
	cout << "========================================" << endl;
	cout << "Total running time: " << (double)(sec_end - sec_start) / CLOCKS_PER_SEC << " [s]" << endl;
	cout << "========================================" << endl;

	print_stop(); // Print time, date and hour when it stops

}

unsigned long long keyFromRunLumiEvent(UInt_t run, UInt_t lumi, ULong64_t event){

  const unsigned long long runMult = 1;
  const unsigned long long lumiMult = 1000000;
  const unsigned long long evtMult = 10000000000;
  const unsigned long long evtLimit = 10000000000;

  unsigned long long key = 0;
  if(event >= evtLimit){
    std::cout << "RUNLUMIEVENTKEY WARNING : \'" << event << "\' is greated that event limit 10^10. returning key 0" << std::endl;
    return key;
  }

  key += runMult* static_cast<unsigned long long>(run);
  key += lumiMult* static_cast<unsigned long long>(lumi);
  key += evtMult*event;

  //std::cout << "key = " << key << std::endl;
  
  return key;
  
}

int main(int argc, char** argv){
				TString firstArgument(argv[1]);
				TString secArgument(argv[2]);				
				TString outfile(argv[3]);
				V0Jet_pPbSkim_ZDC(firstArgument,secArgument,outfile);
}
