#include <TSystem>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void load(){
  //Load all the System libraries
  gSystem->Load("StarRoot");
  gSystem->Load("St_base");
  gSystem->Load("StarClassLibrary"); //problem with vector<double>
  gSystem->Load("StUtilities");

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  
  gSystem->Load("StTpcDb");
  gSystem->Load("StEvent");
  gSystem->Load("StMcEvent");
  gSystem->Load("StMcEventMaker");
  gSystem->Load("StDaqLib");
  gSystem->Load("libgen_Tables");
  gSystem->Load("libsim_Tables");
  gSystem->Load("libglobal_Tables");
  gSystem->Load("StMagF");
  
  gSystem->Load("St_g2t.so");
  gSystem->Load("St_geant_Maker.so");
  gSystem->Load("StAssociationMaker");
  gSystem->Load("StMcAnalysisMaker");
  gSystem->Load("libgeometry_Tables");   
  gSystem->Load("StTriggerUtilities");

  gSystem->Load("StEmcUtil");
  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StPreEclMaker");
  gSystem->Load("StEpcMaker");
  gSystem->Load("StEmcSimulatorMaker");
  
  gSystem->Load("StDbLib");
  gSystem->Load("StDbUtilities");
  gSystem->Load("StDbBroker");
  gSystem->Load("StDetectorDbMaker");
  gSystem->Load("St_db_Maker");
  
  gSystem->Load("StMtdHitMaker");
  gSystem->Load("StMtdUtil");
  gSystem->Load("StMtdMatchMaker");
  gSystem->Load("StMtdCalibMaker");
  gSystem->Load("StBTofUtil");
  gSystem->Load("StVpdCalibMaker");
  
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");

  gSystem->Load("StFmsUtil");
  gSystem->Load("StFmsDbMaker");

  gSystem->Load("StUPCTreeMaker");

  gSystem->Load("StEEmcUtil");
  gSystem->Load("StEEmcDbMaker");
}

//void doEvent(Int_t nEvents=-1, const Char_t *inputFile="root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/AuAu54_production_2017/ReversedFullField/P18ic/2017/170/18170012/st_physics_adc_18170012_raw_1500006.picoDst.root", const TString outputFile="test.root", const Bool_t debug = kFALSE)
//{
void doEvent(Int_t nEvents=-1, const Char_t *inputFile="test.list", const TString outputFile="test.root", const Bool_t debug = kFALSE)
{
	load();

	StChain *chain = new StChain();
	chain->SetDebug(0);
	StMuDebug::setLevel(0); // switch of some debug output

	StMuTimer timer;
	timer.start();

	Bool_t iMuDst = 0;

	string filelist=inputFile;
        size_t npos=filelist.size();
        filelist.insert(npos-5,"new");
        cout<<filelist<<endl;

        char inputFile2[10000];
        sprintf(inputFile2, "%s",filelist.data());
        ofstream fout(inputFile2);

        ifstream fin(inputFile);

        string name;
  
        vector<string> temp_name;
        while(!fin.eof()){
          fin>>name;
			temp_name.push_back( name ); 
        }
        fin.close();
        int number_of_files = temp_name.size();//no bug fix for repeating the last file.	
        for(int i = 0; i < number_of_files; i++){
		fout<<temp_name[i]<<endl;
	}
        fout.close();	
	
	/*	
	ifstream infile;
	infile.open(inputFile);
	string name;
	getline(infile,name);
	infile.close();
        */
	
	std::size_t found = name.find("MuDst.root");
	if(found!=std::string::npos)
		iMuDst = 1;

	if(iMuDst){
		char theFilter[80];
		sprintf(theFilter,".MuDst.root:MuDst.root");
		StMuDstMaker *microMaker = new StMuDstMaker(0,0,"",inputFile2,theFilter,1000);
		microMaker->Init();
		microMaker->SetStatus("*",1);

		
		StEmcADCtoEMaker *adc2e = new StEmcADCtoEMaker();
		adc2e->setPrint(false);
		adc2e->saveAllStEvent(true);//Set to kTRUE if all hits are to be saved on StEvent

		StPreEclMaker *pre_ecl=new StPreEclMaker();
		pre_ecl->setPrint(kFALSE);

		StEpcMaker *epc=new StEpcMaker();
		epc->setPrint(kFALSE);
	}else{
        StPicoDstMaker *picoMaker = new StPicoDstMaker(StPicoDstMaker::IoRead,inputFile2,"picoDst"); 
	}

	St_db_Maker *dbMk = new St_db_Maker("StarDb", "MySQL:StarDb", "$STAR/StarDb","StarDb");
	if(!iMuDst) dbMk->SetDateTime(20160101,0); //for run16 picoDst


	StFmsDbMaker* fmsdb = new StFmsDbMaker("fmsDb");

	StEEmcDbMaker *eemcdb = new StEEmcDbMaker("EEmcDBMaker");

	StUPCTreeMaker *miniTreeMaker = new StUPCTreeMaker();
	miniTreeMaker->setOutFileName(outputFile);
	miniTreeMaker->setFillTree(1);
	//miniTreeMaker->setMaxVtxR(2.);
	//miniTreeMaker->setMaxVtxZ(100.);
	//miniTreeMaker->setMaxVzDiff(3.);
	//miniTreeMaker->setPrintMemory(1);
	//miniTreeMaker->setPrintCpu(1);
	//miniTreeMaker->setPrintConfig(1);
	if(name.find("st_hlt")!=std::string::npos)
		miniTreeMaker->setStreamName("st_hlt");
	if(name.find("st_mtd")!=std::string::npos)
		miniTreeMaker->setStreamName("st_mtd");
	if(name.find("st_ssdmb")!=std::string::npos)
		miniTreeMaker->setStreamName("st_ssdmb");
	if(name.find("st_upc")!=std::string::npos)
		miniTreeMaker->setStreamName("st_upc");
	if(name.find("st_physics")!=std::string::npos)
		miniTreeMaker->setStreamName("st_physics");
	if(debug)
		miniTreeMaker->SetDebug(1);

	if(chain->Init()==kStERR) return;
	cout<<"chain->Init();"<<endl;

	if(nEvents<0){
		if(iMuDst)
			nEvents = microMaker->chain()->GetEntries();
		else 
			nEvents = picoMaker->chain()->GetEntries();
	}

	cout << "****************************************** " << endl;
	cout << "total number of events  " << nEvents << endl;
	cout << "****************************************** " << endl;

	for(Int_t i=0; i<nEvents; i++) {
		if(debug) {
			cout<<endl;
			cout<<"Working on eventNumber "<< i <<endl;
		} else {
			if(i%1000==0)
				cout << "Working on eventNumber " << i << endl;
		}

		chain->Clear();
		int iret = chain->Make(i);

		if(iret) { cout << "Bad return code!" << iret << endl; break;}
	}

	chain->Finish();
	delete chain;

	timer.stop();
	cout << "Total time = " << timer.elapsedTime() << " s" << endl;

	cout << "****************************************** " << endl;
	cout << "Work done... now its time to close up shop!"<< endl;
	cout << "****************************************** " << endl;
}
