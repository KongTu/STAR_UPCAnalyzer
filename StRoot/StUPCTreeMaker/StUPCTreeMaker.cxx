#include "headers.h"
#include "StUPCTreeMaker.h"

ClassImp(StUPCTreeMaker)

//_____________________________________________________________________________
StUPCTreeMaker::StUPCTreeMaker(const Char_t *name) : StMaker(name), 
mFillTree(0), mFillHisto(1), mPrintConfig(1), mPrintMemory(0), mPrintCpu(0), 
mStreamName("st_physics"), fOutFile(0), mOutFileName(""), mEvtTree(0), mVn(2), 
mMaxVtxR(2.0), mMaxVtxZ(30.0), mMaxVzDiff(3.0), mMinTrkPt(0.1), mMaxTrkEta(2.), 
mMinNHitsFit(15), mMinNHitsFitRatio(0.52), mMinNHitsDedx(10), mMaxDca(3.), 
mMaxnSigmaE(2.5), mMaxBeta2TOF(0.03),mEmcCollection(nullptr), mEmcPosition(nullptr), 
mEmcGeom{},mEmcIndex{}
{
  //~
}
//_____________________________________________________________________________
StUPCTreeMaker::~StUPCTreeMaker()
{
	// default destructor
}
//_____________________________________________________________________________
Int_t StUPCTreeMaker::Init()
{
  
  if(!mOutFileName.Length()){
    LOG_ERROR << "StUPCTreeMaker:: no output file specified for tree and histograms." << endm;
    return kStERR;
  }
  fOutFile = new TFile(mOutFileName.Data(),"recreate");
  LOG_INFO << "StUPCTreeMaker:: create the output file to store the tree and histograms: " << mOutFileName.Data() << endm;
  
  if(mFillTree)    bookTree();
  if(mFillHisto)   bookHistos();

  mStPhysics_TriggerIDs.clear();

  if(mStreamName.EqualTo("st_physics")){
    cout<<"add the MB trigger to st_physics"<<endl;

    mStPhysics_TriggerIDs.push_back(580001);
    mStPhysics_TriggerIDs.push_back(580021);
    
  }
  else if(mStreamName.EqualTo("st_ssdmb")){
    cout<<"add the MB trigger to st_ssdmb"<<endl;
    
    mStPhysics_TriggerIDs.push_back(500001);
  }
  
  //initEmc();
  
  return kStOK;
}
//_____________________________________________________________________________
Int_t StUPCTreeMaker::InitRun(const Int_t runnumber)
{
	LOG_INFO<<"Grab runnumber: "<<runnumber<<endm;

	return kStOK;
}
//_____________________________________________________________________________
Int_t StUPCTreeMaker::Finish()
{

	if(fOutFile){
		fOutFile->cd();
		fOutFile->Write();
		fOutFile->Close();
		LOG_INFO << "StUPCTreeMaker::Finish() -> write out tree in " << mOutFileName.Data() << endm;
	}

	if(mPrintConfig) printConfig();

	//finishEmc();
	return kStOK;
}
//_____________________________________________________________________________
Int_t StUPCTreeMaker::Make()
{
	memset(&mEvtData, 0, sizeof(mEvtData)); //initial the mEvtData structure

	StTimer timer;
	if(mPrintMemory) StMemoryInfo::instance()->snapshot();
	if(mPrintCpu)    timer.start();

	mMuDstMaker = (StMuDstMaker *)GetMaker("MuDst");
	mPicoDstMaker = (StPicoDstMaker *)GetMaker("picoDst");
  fmsDbMaker = (StFmsDbMaker *)GetMaker("fmsDb");

	if(Debug()){
		LOG_INFO<<"MuDstMaker pointer: "<<mMuDstMaker<<endm;
		LOG_INFO<<"PicoDstMaker pointer: "<<mPicoDstMaker<<endm;
	}

	if(mMuDstMaker){
		if(Debug()) LOG_INFO<<"Use MuDst file as input"<<endm;
		mMuDst = mMuDstMaker->muDst();
		if(!mMuDst){
			LOG_WARN<<"No MuDst !"<<endm;
			return kStOK;
		}
	}
	else if(mPicoDstMaker){
		if(Debug()) LOG_INFO<<"Use Pico file as input"<<endm;
		mPicoDst = mPicoDstMaker->picoDst();
		if(!mPicoDst){
			LOG_WARN<<"No PicoDst !"<<endm;
			return kStOK;
		}
	}
	else{
		LOG_WARN<<"No StMuDstMaker and No StPicoDstMaker !"<<endm;
		return kStOK;
	}

	if(mPicoDstMaker){
        if(!processPicoEvent())  return kStOK;
	}

	if(mFillTree) mEvtTree->Fill();

	if(mPrintMemory){
		StMemoryInfo::instance()->snapshot();
		StMemoryInfo::instance()->print();
	}

	if(mPrintCpu){
		timer.stop();
		LOG_INFO << "CPU time for StUPCTreeMaker::Make(): " 
			<< timer.elapsedTime() << "sec " << endm;
	}

	return kStOK;
}
//_____________________________________________________________________________
Bool_t StUPCTreeMaker::processPicoEvent()
{
  double PCAPtBin[10] = {0.2,0.5,0.8,1.2,1.6,2,2.5,3,4,6};
  // double etaBin[] = {-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};

  if(mFillHisto) hEvent->Fill(0.5);
  
  StPicoEvent *picoEvent = mPicoDst->event();
  if(!picoEvent){
    LOG_WARN<<"No event level information !"<<endm;
    return kFALSE;
  }
  
  Bool_t validTrigger = kFALSE;
  Bool_t VPD5 = kFALSE;
  Bool_t VPD5HM = kFALSE;

  Int_t nTrigs = 0;
  if(mStreamName.EqualTo("st_physics")){
    for(Int_t i=0;i<mStPhysics_TriggerIDs.size();i++){
      //cout<<"******** "<<mStPhysics_TriggerIDs[0]<<endl;
      if(picoEvent->isTrigger(mStPhysics_TriggerIDs[i])){
	validTrigger = kTRUE;
	
        mTrigId[nTrigs] = mStPhysics_TriggerIDs[i];

	if(i==0) VPD5 = kTRUE;
	if(i==1) VPD5HM = kTRUE;

	if(mFillHisto) hEvent->Fill(0.5+i+1);
	nTrigs++;
      }
    }
  }
  else{
    LOG_WARN<<"The data stream name is wrong !"<<endm;
    return kFALSE;
  }
  
  mNTrigs = nTrigs;
  
  if(!validTrigger){
    //LOG_WARN<<"No valid MB related triggers !"<<endm;
    return kFALSE;
  }
  
  mRunId          = picoEvent->runId();
  mEventId        = picoEvent->eventId();
  mRefMult        = picoEvent->refMult();
  mGRefMult       = picoEvent->grefMult();
  mBBCRate        = picoEvent->BBCx();
  mZDCRate        = picoEvent->ZDCx();
  mBField         = picoEvent->bField();
  mVpdVz          = picoEvent->vzVpd();
  
  StThreeVectorF vtxPos    = picoEvent->primaryVertex();
  mVertexX        = vtxPos.x();
  mVertexY        = vtxPos.y();
  mVertexZ        = vtxPos.z();
  if(Debug()){
    LOG_INFO<<"RunId: "<<mRunId<<endm;
    LOG_INFO<<"EventId: "<<mEventId<<endm;
    LOG_INFO<<"mZDCX: "<<mZDCRate<<endm;
    LOG_INFO<<"VPD Vz: "<<mVpdVz<<" \tTPC Vz: "<<mVertexZ<<endm;
  }

  if(mFillHisto){
    hVtxYvsVtxX->Fill(mVertexX, mVertexY);
    hVPDVzvsTPCVz->Fill(mVertexZ, mVpdVz);
    hVzDiff->Fill(mVertexZ - mVpdVz);
  }

  if(TMath::Abs(vtxPos.x())<1.e-5 && TMath::Abs(vtxPos.y())<1.e-5 && TMath::Abs(vtxPos.z())<1.e-5) return kFALSE;
  if(mFillHisto) hEvent->Fill(3.5);
  if(sqrt(vtxPos.x()*vtxPos.x()+vtxPos.y()*vtxPos.y())>=mMaxVtxR) return kFALSE;
  if(mFillHisto) hEvent->Fill(4.5);
  if(TMath::Abs(vtxPos.z())>=mMaxVtxZ) return kFALSE;
  if(mFillHisto) hEvent->Fill(5.5);
  if(TMath::Abs(mVertexZ - mVpdVz)>=mMaxVzDiff) return kFALSE;
  if(mFillHisto) hEvent->Fill(6.5);
  
  if(mFillHisto){
    if(VPD5) hEvent->Fill(7.5);
    if(VPD5HM) hEvent->Fill(8.5);
  }
  
  hVtxZ->Fill( vtxPos.z() );


  mNTofHits = mPicoDst->numberOfBTofHits();//better one for multiplicity
  mZDCeast = picoEvent->ZdcSumAdcEast();
  mZDCwest = picoEvent->ZdcSumAdcWest();
  
  for(int ch=0;ch<24;ch++) {
    mBbcQ[ch]    = picoEvent->bbcAdcEast(ch);
    mBbcQ[ch+24] = picoEvent->bbcAdcWest(ch);

  }

  Int_t nNodes = mPicoDst->numberOfTracks();
  if(Debug()){
    LOG_INFO<<"# of global tracks in picoDst: "<<nNodes<<endm;
  }

  memset(mNHitsFit, 0, sizeof(mNHitsFit));
  memset(mNHitsPoss, 0, sizeof(mNHitsPoss));
  //centrality
  Int_t   NCentbin = 9;
  Int_t   CentralityBins  [9] = { 377, 327, 248, 186, 138, 99, 68, 45, 26} ;
  Int_t   MiddleBinID     [9] = { 0,1,2,3,4,5,6,7,8};  // ID Number
    
  Int_t   myCentrality=-1;

  hRefMult->Fill( mRefMult );
  
  for(int i=0;i!=9;i++){
        if( mRefMult > CentralityBins[i] ){
            myCentrality = MiddleBinID[i];
            break;
        }
        else{
            myCentrality = -1;
        }
    }
    
  if(myCentrality<0) return kFALSE;

  hCentrality->Fill( myCentrality );

  Int_t nTrks    = 0;

  //Q-vectors
  TComplex Q_n1_pt[20][2], Q_0_pt[20][2];//pt 1 dimention
  TComplex Q_n3_1_FMSplus, Q_0_1_FMSplus;
  TComplex Q_n3_1_TPCminus, Q_0_1_TPCminus;
  TComplex Q_n3_1_TPCplus, Q_0_1_TPCplus;

  //loop FMS
  unsigned int nFms = mPicoDst->numberOfFmsHits();
  for(unsigned int ihit=0; ihit<nFms; ihit++){
      StPicoFmsHit *pFms = mPicoDst->fmsHit(ihit);
      
      if(pFms->adc()>0){
          int id = -999;
          if(pFms->detectorId()==10||pFms->detectorId()==11) id=1;//inner
          if(pFms->detectorId()==8 ||pFms->detectorId()==9 ) id=0;//outter
          
          if(id<0) continue;
          
          short c=pFms->channel();
          short d=pFms->detectorId();
          
          StThreeVectorF xyz = fmsDbMaker->getStarXYZ(d,c);//pFms->starXYZ();
          
          float x = xyz.x();
          float y = xyz.y();
          float z = xyz.z();
          
          if(z<600) continue;
          
          float r = sqrt(x*x + y*y);
          float phi = atan2(y, x);
          
          float the0 = atan2(r,z - mVertexZ);
          
          hFmsXYdis->Fill(x,y);
          
          float energy = fmsDbMaker->getGain(d,c)*fmsDbMaker->getGainCorrection(d,c)*pFms->adc();
          
          double Et = energy*sin(the0);
          
          if(Et<0.01) continue;

          double w = Et;
          Q_n3_1_FMSplus += q_vector(mVn, 1, w, phi);
          Q_0_1_FMSplus += q_vector(0, 1, w, phi);
          
      }
  }//end of FMS

  //track loop:
  for(Int_t i=0;i<nNodes;i++){
    
    StPicoTrack *pTrack = mPicoDst->track(i);
    if(!pTrack) continue;
        
    //track cut
    if(!isValidTrack(pTrack, vtxPos)) continue; 
    
    mPt[nTrks] = -999;
    mEta[nTrks] = -999;
    mPhi[nTrks] = -999;
    
    mgPt[nTrks] = -999;
    mgEta[nTrks] = -999;
    mgPhi[nTrks] = -999;

    mCharge[nTrks] = -999;

    mCharge[nTrks] = pTrack->charge();
    
    StThreeVectorF pMom = pTrack->pMom();
    StThreeVectorF gMom = pTrack->gMom();
    StThreeVectorF origin = pTrack->origin();
    
    mPt[nTrks] = pMom.perp();
    mEta[nTrks] = pMom.pseudoRapidity();
    mPhi[nTrks] = pMom.phi();
    
    mgPt[nTrks]              = gMom.perp();
    mgEta[nTrks]             = gMom.pseudoRapidity();
    mgPhi[nTrks]             = gMom.phi();
    //		mgOriginX[nTrks]         = origin.x();
    //		mgOriginY[nTrks]         = origin.y();
    //		mgOriginZ[nTrks]         = origin.z();

    mNHitsFit[nTrks]         = pTrack->nHitsFit();
    mNHitsPoss[nTrks]        = pTrack->nHitsMax();
    //mNHitsDedx[nTrks]        = pTrack->nHitsDedx();
    mDedx[nTrks]             = pTrack->dEdx(); 
    //mDndx[nTrks]             = pTrack->dNdx();
    //mDndxError[nTrks]        = pTrack->dNdxError();
    mNSigmaE[nTrks]          = pTrack->nSigmaElectron();
    //mDca[nTrks]              = (pTrack->dca()-vtxPos).mag();
    mDca[nTrks]              = (pTrack->dcaPoint()-vtxPos).mag();
    //mIsHFTTrk[nTrks]         = pTrack->isHFTTrack();
    //mHasHFT4Layers[nTrks]    = pTrack->hasHft4Layers();
    
    if(mFillHisto){
    hdEdxvsP->Fill(pMom.mag(), mDedx[nTrks]);
    hdNdxvsP->Fill(pMom.mag(), mDndx[nTrks]);
    hnSigEvsP->Fill(pMom.mag(), mNSigmaE[nTrks]);
    }
    
    Int_t bTofPidTraitsIndex = pTrack->bTofPidTraitsIndex();
    //mTOFMatchFlag[nTrks]     = -1; 
    mTOFLocalY[nTrks]        = -999;
    mBeta2TOF[nTrks]         = -999;
    if(bTofPidTraitsIndex>=0){
      StPicoBTofPidTraits *btofPidTraits = mPicoDst->btofPidTraits(bTofPidTraitsIndex);
      //mTOFMatchFlag[nTrks] = btofPidTraits->btofMatchFlag(); 
      mTOFLocalY[nTrks]    = btofPidTraits->btofYLocal();
      mBeta2TOF[nTrks]     = btofPidTraits->btofBeta();
      
      if(mFillHisto) hBetavsP->Fill(pMom.mag(), 1./mBeta2TOF[nTrks]);
    }
    
    //add emc matching

    Int_t bemcPidTraitsIndex              = pTrack->bemcPidTraitsIndex();
    if( bemcPidTraitsIndex>=0 ){
      StPicoBEmcPidTraits *bemcPidTraits = mPicoDst->bemcPidTraits(bemcPidTraitsIndex);
      mBEMCE[nTrks]         = bemcPidTraits->bemcE();
    }

    mNEmc=0;

    nTrks++;

    //calculate Qvectors
    double pt  = pTrack->pMom().perp();
    double trkEta = pTrack->pMom().pseudoRapidity();
    double phi = pTrack->pMom().phi();
    double dca = (pTrack->dcaPoint()-vtxPos).mag();
    
    if(TMath::Abs(trkEta)>1.) continue;
    if(pTrack->nHitsFit()<15) continue;
    if(pTrack->nHitsFit()*1./pTrack->nHitsMax()<0.52) continue;//zhenyu's version is 0.52
    if(dca>3.) continue;
    
    double weight = 1.0;

    if( pt > 0.4 && pt < 2.0 ){
      //QnB and QnC event plane Q vectors:
      if( trkEta < -0.5 && trkEta > -1.0 ){

        Q_n3_1_TPCminus += q_vector(mVn, 1, weight, phi);
        Q_0_1_TPCminus += q_vector(0, 1, weight, phi);

      }
      if( trkEta > 0.5 && trkEta < 1.0 ){

        Q_n3_1_TPCplus += q_vector(mVn, 1, weight, phi);
        Q_0_1_TPCplus += q_vector(0, 1, weight, phi);

      }
    }

    if(pt < 0.2) continue;//for POI pt selections

    /*start pT*/
    for(int ipt = 0; ipt < 9; ipt++){
      if( pt > PCAPtBin[ipt] && pt < PCAPtBin[ipt+1] ){

        if( pTrack->charge() == +1 ){//positive charge

            Q_n1_pt[ipt][0] += q_vector(mVn, 1, weight, phi);
            Q_0_pt[ipt][0] += q_vector(0, 1, weight, phi);
        }
        if( pTrack->charge() == -1 ){//negative charge

            Q_n1_pt[ipt][1] += q_vector(mVn, 1, weight, phi);
            Q_0_pt[ipt][1] += q_vector(0, 1, weight, phi);
        }

      }//end if
    }//end pT loop
    //end pT

  }//end of track loop
  

  //put centrality cut here. 

/*Three sub-events resolution*/

  for(int icent = 0; icent < NCentbin; icent++){

    if(myCentrality!=icent) continue;

    TComplex N_2_trk, D_2_trk;
    N_2_trk = Q_n3_1_TPCplus*TComplex::Conjugate(Q_n3_1_TPCminus);
    D_2_trk = Q_0_1_TPCplus*Q_0_1_TPCminus;

    cn_QbQc[icent]->Fill(N_2_trk.Re()/D_2_trk.Re(), D_2_trk.Re());

    N_2_trk = Q_n3_1_FMSplus*TComplex::Conjugate(Q_n3_1_TPCminus);
    D_2_trk = Q_0_1_FMSplus*Q_0_1_TPCminus;

    cn_QaQb[icent]->Fill(N_2_trk.Re()/D_2_trk.Re(), D_2_trk.Re());

    N_2_trk = Q_n3_1_FMSplus*TComplex::Conjugate(Q_n3_1_TPCplus);
    D_2_trk = Q_0_1_FMSplus*Q_0_1_TPCplus;

    cn_QaQc[icent]->Fill(N_2_trk.Re()/D_2_trk.Re(), D_2_trk.Re());

    cn_Qa_real[icent]->Fill(Q_n3_1_FMSplus.Re()/Q_0_1_FMSplus.Re(), Q_0_1_FMSplus.Re() );
    cn_Qa_imag[icent]->Fill(Q_n3_1_FMSplus.Im()/Q_0_1_FMSplus.Re(), Q_0_1_FMSplus.Re() );
   
    cn_Qb_real[icent]->Fill(Q_n3_1_TPCminus.Re()/Q_0_1_TPCminus.Re(), Q_0_1_TPCminus.Re() );
    cn_Qb_imag[icent]->Fill(Q_n3_1_TPCminus.Im()/Q_0_1_TPCminus.Re(), Q_0_1_TPCminus.Re() );

    cn_Qc_real[icent]->Fill(Q_n3_1_TPCplus.Re()/Q_0_1_TPCplus.Re(), Q_0_1_TPCplus.Re() );
    cn_Qc_imag[icent]->Fill(Q_n3_1_TPCplus.Im()/Q_0_1_TPCplus.Re(), Q_0_1_TPCplus.Re() );


  /*pt vn*/

    for(int ipt = 0; ipt < 9; ipt++){
      for(int charge = 0; charge < 2; charge++){

        TComplex N_2_trk, D_2_trk;
        N_2_trk = Q_n1_pt[ipt][charge] * TComplex::Conjugate( Q_n3_1_FMSplus );
        D_2_trk = Q_0_pt[ipt][charge]*Q_0_1_FMSplus;

        cn_tracker_fms[icent][ipt][charge]->Fill(N_2_trk.Re()/D_2_trk.Re(), D_2_trk.Re());

        cn_tracker_fms_real[icent][ipt][charge]->Fill(Q_n1_pt[ipt][charge].Re()/Q_0_pt[ipt][charge].Re(), Q_0_pt[ipt][charge].Re());
        cn_tracker_fms_imag[icent][ipt][charge]->Fill(Q_n1_pt[ipt][charge].Im()/Q_0_pt[ipt][charge].Re(), Q_0_pt[ipt][charge].Re());

      }
    }
  } 


  mNTrks       = nTrks;

  if(Debug()){
    LOG_INFO<<"# of primary tracks stored: "<<mNTrks<<endm;
    //LOG_INFO<<"# of EMC matched Tracks stored: "<<mNBEMCTrks<<endm;
    //LOG_INFO<<"# of MTD matched Tracks stored: "<<mNMTDTrks<<endm;
  }
    
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t StUPCTreeMaker::isValidTrack(StPicoTrack *pTrack, StThreeVectorF vtxPos) const
{
	//Float_t pt  = pTrack->pMom().perp();
	//Float_t eta = pTrack->pMom().pseudoRapidity();
	//Float_t dca = (pTrack->dca()-vtxPos).mag();
	Float_t dca = (pTrack->dcaPoint()-vtxPos).mag();

	//if(pt<mMinTrkPt)                            return kFALSE;
	//if(TMath::Abs(eta)>mMaxTrkEta)              return kFALSE;
	if(pTrack->nHitsFit()<mMinNHitsFit)         return kFALSE;
	if(pTrack->nHitsFit()*1./pTrack->nHitsMax()<mMinNHitsFitRatio) return kFALSE;
	//if(pTrack->nHitsDedx()<mMinNHitsDedx)       return kFALSE;
	if(dca>mMaxDca)                             return kFALSE;

	return kTRUE;
}
TComplex StUPCTreeMaker::q_vector(double n, double p, double w, double phi) 
{
  double term1 = pow(w,p);
  TComplex e(1, n*phi, 1);
  return term1*e;
}
//_____________________________________________________________________________
void StUPCTreeMaker::bookTree()
{
	LOG_INFO << "StUPCTreeMaker:: book the event tree to be filled." << endm;

	mEvtTree = new TTree("miniDst","miniDst for smallsystem");
	mEvtTree->SetAutoSave(100000000); // 100 MB

	// event information
	mEvtTree->Branch("mRunId", &mRunId, "mRunId/I");
	mEvtTree->Branch("mEventId", &mEventId, "mEventId/I");
	
	mEvtTree->Branch("mNTrigs", &mNTrigs, "mNTrigs/I");
	mEvtTree->Branch("mTrigId", mTrigId, "mTrigId[mNTrigs]/I");
	
	mEvtTree->Branch("mRefMult", &mRefMult, "mRefMult/S");
	mEvtTree->Branch("mGRefMult", &mGRefMult, "mGRefMult/S");
	
	mEvtTree->Branch("mBBCRate", &mBBCRate, "mBBCRate/F");
	mEvtTree->Branch("mZDCRate", &mZDCRate, "mZDCRate/F");
	mEvtTree->Branch("mBField", &mBField, "mBField/F");
	mEvtTree->Branch("mVpdVz", &mVpdVz, "mVpdVz/F");
	mEvtTree->Branch("mVertexX", &mVertexX, "mVertexX/F");
	mEvtTree->Branch("mVertexY", &mVertexY, "mVertexY/F");
	mEvtTree->Branch("mVertexZ", &mVertexZ, "mVertexZ/F");


	//BBC
	mEvtTree->Branch("mBbcQ", mBbcQ, "mBbcQ[48]/I");

	//ZDC
	mEvtTree->Branch("mZDCwest", &mZDCwest, "mZDCwest/I");
	mEvtTree->Branch("mZDCeast", &mZDCeast, "mZDCeast/I");

	//nTOF
  mEvtTree->Branch("mNTofHits",&mNTofHits,"mNTofHits/I");

	//Emc
  //mEvtTree->Branch("mNEmc",    &mNEmc,     "mNEmc/I");
	//mEvtTree->Branch("mEmcE",     mEmcE,     "mEmcE[mNEmc]/F");
	//mEvtTree->Branch("mEmcEta",   mEmcEta,   "mEmcEta[mNEmc]/F");
	//mEvtTree->Branch("mEmcPhi",   mEmcPhi,   "mEmcPhi[mNEmc]/F");

	//all tracks information
	
	mEvtTree->Branch("mNTrks", &mNTrks, "mNTrks/I");
	//mEvtTree->Branch("mTrkId", mTrkId, "mTrkId[mNTrks]/S");

	mEvtTree->Branch("mCharge", mCharge, "mCharge[mNTrks]/I");
	mEvtTree->Branch("mPt", mPt, "mPt[mNTrks]/F");
	mEvtTree->Branch("mEta", mEta, "mEta[mNTrks]/F");
	mEvtTree->Branch("mPhi", mPhi, "mPhi[mNTrks]/F");

	mEvtTree->Branch("mgPt", mgPt, "mgPt[mNTrks]/F");
  mEvtTree->Branch("mgEta", mgEta, "mgEta[mNTrks]/F");
  mEvtTree->Branch("mgPhi", mgPhi, "mgPhi[mNTrks]/F");

	mEvtTree->Branch("mNHitsFit", mNHitsFit, "mNHitsFit[mNTrks]/I");
	mEvtTree->Branch("mNHitsPoss", mNHitsPoss, "mNHitsPoss[mNTrks]/I");

	mEvtTree->Branch("mDca", mDca, "mDca[mNTrks]/F");

	mEvtTree->Branch("mTOFLocalY", mTOFLocalY, "mTOFLocalY[mNTrks]/F");
	mEvtTree->Branch("mBeta2TOF", mBeta2TOF, "mBeta2TOF[mNTrks]/F");
	
	mEvtTree->Branch("mBEMCE", mBEMCE, "mBEMCE[mNTrks]/F");
	//mEvtTree->Branch("mBEMCE0", mBEMCE0, "mBEMCE0[mNTrks]/F");
	//mEvtTree->Branch("mBEMCE1", mBEMCE1, "mBEMCE1[mNTrks]/F");
}
//_____________________________________________________________________________
void StUPCTreeMaker::bookHistos()
{
    //shengli's event keeping
	hEvent = new TH1D("hEvent","Event statistics",25,0,25);
	hEvent->GetXaxis()->SetBinLabel(1,"All events");
	hEvent->GetXaxis()->SetBinLabel(2,"VPD5");
	hEvent->GetXaxis()->SetBinLabel(3,"VPD5HM");
	hEvent->GetXaxis()->SetBinLabel(4,"None-Zero Vertex");
	hEvent->GetXaxis()->SetBinLabel(5,Form("|V_{r}|<%1.2f cm",mMaxVtxR));
	hEvent->GetXaxis()->SetBinLabel(6,Form("|V_{z}|<%1.2f cm",mMaxVtxZ));
	hEvent->GetXaxis()->SetBinLabel(7,Form("|V_{z}Diff|<%1.2f cm",mMaxVzDiff));
	hEvent->GetXaxis()->SetBinLabel(8,"VPD5");
	hEvent->GetXaxis()->SetBinLabel(9,"VPD5HM");

  hVtxZ = new TH1D("hVtxZ","hVtxZ",1000,-50,50);
  hRefMult = new TH1D("hRefMult","hRefMult",50000,0,50000);

	hVtxYvsVtxX = new TH2D("hVtxYvsVtxX","hVtxYvsVtxX; V_{x} (cm); V_{y} (cm)",120,-3,3,120,-3,3); 
	hVPDVzvsTPCVz = new TH2D("hVPDVzvsTPCVz","hVPDVzvsTPCVz; TPC V_{z} (cm); VPD V_{z} (cm)",200,-50,50,200,-50,50);
	hVzDiff = new TH1D("hVzDiff","hVzDiff; TPC V_{z} - VPD V_{z} (cm)",80,-20,20);
	hGRefMultvsGRefMultCorr = new TH2D("hGRefMultvsGRefMultCorr","hGRefMultvsGRefMultCorr; grefMultCorr; grefMult",1000,0,1000,1000,0,1000);
	hCentrality = new TH1D("hCentrality","hCentrality; mCentrality",16,0,16);
    
  //Q vector histograms:
  for(int icent = 0; icent < 9; icent++){

    cn_Qa_real[icent] = new TH1D(Form("cn_Qa_real_%d",icent),"cn_Qa_real",100,-1,1);
    cn_Qb_real[icent] = new TH1D(Form("cn_Qb_real_%d",icent),"cn_Qb_real",100,-1,1);
    cn_Qc_real[icent] = new TH1D(Form("cn_Qc_real_%d",icent),"cn_Qc_real",100,-1,1);
    
    cn_Qa_imag[icent] = new TH1D(Form("cn_Qa_imag_%d",icent),"cn_Qa_imag",100,-1,1);
    cn_Qb_imag[icent] = new TH1D(Form("cn_Qb_imag_%d",icent),"cn_Qb_imag",100,-1,1);
    cn_Qc_imag[icent] = new TH1D(Form("cn_Qc_imag_%d",icent),"cn_Qc_imag",100,-1,1);

    cn_QaQb[icent] = new TH1D(Form("cn_QaQb_%d",icent),"cn_QaQb",100,-1,1);
    cn_QaQc[icent] = new TH1D(Form("cn_QaQc_%d",icent),"cn_QaQc",100,-1,1);
    cn_QbQc[icent] = new TH1D(Form("cn_QbQc_%d",icent),"cn_QbQc",100,-1,1);

    for(int ipt = 0; ipt < 9; ipt++){
      for(int charge = 0; charge < 2; charge++){

        //cn_tracker[ipt][charge] = new TH1D(Form("cn_tracker_%d_%d", ipt, charge),"cn_tracker",100,-1,1);
        cn_tracker_fms[icent][ipt][charge] = new TH1D(Form("cn_tracker_fms_%d_%d_%d",icent, ipt, charge),"cn_tracker_fms",100,-1,1);
        cn_tracker_fms_real[icent][ipt][charge] = new TH1D(Form("cn_tracker_fms_real_%d_%d_%d",icent, ipt, charge),"cn_tracker_fms",100,-1,1);
        cn_tracker_fms_imag[icent][ipt][charge] = new TH1D(Form("cn_tracker_fms_imag_%d_%d_%d",icent, ipt, charge),"cn_tracker_fms",100,-1,1);

      }
    }
  } 
    
	hdEdxvsP = new TH2D("hdEdxvsP","hdEdxvsP; p (GeV/c); dE/dx (KeV/cm)",300,0,15,400,0,20);
	hdNdxvsP = new TH2D("hdNdxvsP","hdNdxvsP; p (GeV/c); dN/dx",300,0,15,400,0,200);
	hnSigEvsP = new TH2D("hnSigEvsP","hnSigEvsP; p (GeV/c); n#sigma_{e}",300,0,15,700,-15,20);
	hBetavsP = new TH2D("hBetavsP","hBetavsP; p (GeV/c); 1/#beta",300,0,15,800,0,4);
  hFmsXYdis = new TH2D("hFmsXYdis", "hFmsXYdis", 200, -200.0, 200.0, 200, -200, 200);
}
//_____________________________________________________________________________
void StUPCTreeMaker::printConfig()
{
	const char *decision[2] = {"no","yes"};
	printf("=== Configuration for StUPCTreeMaker ===\n");
	printf("Fill the miniDst tree: %s\n",decision[mFillTree]);
	printf("Fill the QA histo: %s\n",decision[mFillHisto]);
	printf("Maximum |Vr|: %1.2f\n",mMaxVtxR);
	printf("Maximum |Vz|: %1.2f\n",mMaxVtxZ);
	printf("Maximum |VzDiff|: %1.2f\n",mMaxVzDiff);
	printf("Minimum Track pt: %1.2f\n",mMinTrkPt);
	printf("Maximum Track |eta| : %1.2f\n",mMaxTrkEta);
	printf("Minimum number of fit hits: %d\n",mMinNHitsFit);
	printf("Minimum ratio of fit hits: %1.2f\n",mMinNHitsFitRatio);
	printf("Minimum number of dedx hits: %d\n",mMinNHitsDedx);
	printf("Maximum dca: %1.2f\n",mMaxDca);
	printf("Maximum |nSigmaE| for TPCe: %1.2f\n",mMaxnSigmaE);
	printf("Maximum |1-1/beta| for TPCe: %1.2f\n",mMaxBeta2TOF);
	printf("=======================================\n");
}

//____________________________________________________________________________
/*void StUPCTreeMaker::initEmc()
{
  mEmcPosition = new StEmcPosition();

  for (int i = 0; i < 4; ++i)
  {
    mEmcGeom[i] = StEmcGeom::getEmcGeom(detname[i].Data());
  }
}

void StUPCTreeMaker::finishEmc()
{
  delete mEmcPosition;
  mEmcPosition = nullptr;

  std::fill_n(mEmcGeom, 4, nullptr);
}*/

//_________________
/*void StUPCTreeMaker::buildEmcIndex()
{
  StEmcDetector* mEmcDet = mMuDst->emcCollection()->detector(kBarrelEmcTowerId);
  std::fill_n(mEmcIndex, sizeof(mEmcIndex) / sizeof(mEmcIndex[0]), nullptr);

  if (!mEmcDet) return;
  for (size_t iMod = 1; iMod <= mEmcDet->numberOfModules(); ++iMod)
  {
    StSPtrVecEmcRawHit& modHits = mEmcDet->module(iMod)->hits();
    for (size_t iHit = 0; iHit < modHits.size(); ++iHit)
    {
      StEmcRawHit* rawHit = modHits[iHit];
      if (!rawHit) continue;
      unsigned int softId = rawHit->softId(1);
      if (mEmcGeom[0]->checkId(softId) == 0) // OK
      {
        mEmcIndex[softId - 1] = rawHit;
      }

    }
  }
}*/
