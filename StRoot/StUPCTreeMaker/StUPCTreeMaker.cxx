#include "headers.h"
#include "StUPCTreeMaker.h"

ClassImp(StUPCTreeMaker)

//_____________________________________________________________________________
StUPCTreeMaker::StUPCTreeMaker(const Char_t *name) : StMaker(name), 
mFillTree(0), mFillHisto(1), mPrintConfig(1), mPrintMemory(0), mPrintCpu(0),mDoMC_(1), 
mStreamName("st_upc"), fOutFile(0), mOutFileName(""), mEvtTree(0),
mMaxVtxR(999.0), mMaxVtxZ(130.0), mMinTrkPt(0.2), mMaxTrkEta(1.1), 
mMinNHitsFit(13), mMinNHitsFitRatio(0.50), mMinNHitsDedx(10), mMaxDca(10.), 
mMaxnSigmaE(999), mMaxBeta2TOF(0.99),mEmcCollection(nullptr), mEmcPosition(nullptr), 
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
  mEmcPosition = new StEmcPosition();
  for(Int_t i=0;i<4;i++){
    if(i==1) continue;
    mEmcGeom[i] = StEmcGeom::getEmcGeom(detname[i].Data());
  }

  if(!mOutFileName.Length()){
    LOG_ERROR << "StUPCTreeMaker:: no output file specified for tree and histograms." << endm;
    return kStERR;
  }
  fOutFile = new TFile(mOutFileName.Data(),"recreate");
  LOG_INFO << "StUPCTreeMaker:: create the output file to store the tree and histograms: " << mOutFileName.Data() << endm;
  
  if(mFillTree)    bookTree();
  if(mFillHisto)   bookHistos();

  mStUPC_TriggerIDs.clear();

  if(mStreamName.EqualTo("st_upc")){
    
    cout<<"add the dAu 200 GeV UPC trigger to st_upc"<<endl;
    mStUPC_TriggerIDs.push_back(530701);
    mStUPC_TriggerIDs.push_back(530702);
    mStUPC_TriggerIDs.push_back(530703);
    
  }
  else if(mStreamName.EqualTo("st_ssdmb")){
    cout<<"add the MB trigger to st_ssdmb"<<endl;
    
    mStPhysics_TriggerIDs.push_back(500001);
  }
  
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
    mBemcTriggerSimu = 0;
    mTriggerSimuMaker = (StTriggerSimuMaker*)GetMaker("StarTrigSimu");
    if(mTriggerSimuMaker){
            mBemcTriggerSimu  = (StBemcTriggerSimu*)mTriggerSimuMaker->bemc;
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
  if(mMuDstMaker){
    if(!processMuDstEvent()) return kStOK;
  }
	else{
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
//-----------------------------------------------------------------------------
Bool_t StUPCTreeMaker::processMuDstEvent()
{
  if(mFillHisto) hEvent->Fill(0.5);

  StMuEvent *mMuEvent = mMuDst->event();
  if(!mMuEvent){
    LOG_WARN<<"No event level information !"<<endm;
    return kFALSE;
  }

  Bool_t validTrigger   = kFALSE;

  Int_t nTrigs = 0;
  if(mStreamName.EqualTo("st_upc")){
    for(unsigned i=0;i<mStUPC_TriggerIDs.size();i++){
      if(mMuEvent->triggerIdCollection().nominal().isTrigger(mStUPC_TriggerIDs[i])){
        
        validTrigger   = kTRUE;
        mTrigId[nTrigs] = mStUPC_TriggerIDs[i];
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
    if( !mDoMC_ ){
      LOG_WARN<<"No valid UPC related triggers !"<<endm;
      return kFALSE;
    }
  }

  if(mFillHisto){
    if(validTrigger) {
      hEvent->Fill(6.5);
    }
  }

  mRunId          = mMuEvent->runId();
  mEventId        = mMuEvent->eventId();
  mRefMult        = mMuEvent->refMult();
  mGRefMult       = mMuEvent->grefmult();
  mBBCRate        = mMuEvent->runInfo().bbcCoincidenceRate();
  mZDCRate        = mMuEvent->runInfo().zdcCoincidenceRate();
  mBField         = mMuEvent->runInfo().magneticField();

  if(Debug()){
    LOG_INFO<<"RunId: "<<mRunId<<endm;
    LOG_INFO<<"EventId: "<<mEventId<<endm;
    LOG_INFO<<"mZDCX: "<<mZDCRate<<endm;
    LOG_INFO<<"VPD Vz: "<<mVpdVz<<" \tTPC Vz: "<<mVertexZ<<endm;
  }

  int Nvertex = mMuDst->primaryVertices()->GetEntriesFast();
  hNvertex->Fill( Nvertex );
  int bestvertex = -1;

  for( int jvtx = 0; jvtx < Nvertex; jvtx++){

    StMuDst::setVertexIndex(jvtx);
    //get array of primary tracks
    TObjArray *trkArray = mMuDst->primaryTracks();
    if( !trkArray ) continue;

    int nElectrons = 0;
    int nMatchTof = 0;
    int nTracks = 0;
    Short_t nSelTracks    = 0;
    Short_t nSelBEMCTrks = 0;

    for(Int_t itrk=0; itrk<trkArray->GetEntriesFast(); itrk++) {
      
      StMuTrack *track = dynamic_cast<StMuTrack*>( trkArray->At(itrk) );
      if( !track ) continue;

      //no cuts:
      nTracks++;
      //matching to BEMC cluster
      bool matchBemc = false;
      matchBemc = getBemcInfo(track,nSelTracks,nSelBEMCTrks,true);
    
      //TOF matching
      const StMuBTofPidTraits &tofPid = track->btofPidTraits();
      Bool_t matchTof = tofPid.matchFlag() != 0 ? kTRUE : kFALSE;
      if( matchTof ) nMatchTof++;

      //require at least one match, only in data
      if( matchBemc ) nElectrons++;      
    }

    if( nElectrons >= 2 && nTracks >= 2 ){ bestvertex = jvtx; break;}
    if( nMatchTof >= 2 && nTracks >= 2 ) { bestvertex = jvtx; break;}
 
  }

  hbestVertex->Fill( bestvertex );

  if( bestvertex != -1 ) {StMuDst::setVertexIndex(bestvertex);}

  StThreeVectorF vtxPos    = mMuEvent->primaryVertexPosition();
  mVertexX        = vtxPos.x();
  mVertexY        = vtxPos.y();
  mVertexZ        = vtxPos.z();

  if(TMath::Abs(vtxPos.x())<1.e-5 && TMath::Abs(vtxPos.y())<1.e-5 && TMath::Abs(vtxPos.z())<1.e-5) return kFALSE;
  if(mFillHisto) hEvent->Fill(10.5);
  if(sqrt(vtxPos.x()*vtxPos.x()+vtxPos.y()*vtxPos.y())>=mMaxVtxR) return kFALSE;
  if(mFillHisto) hEvent->Fill(11.5);
  if(TMath::Abs(vtxPos.z())>=mMaxVtxZ) return kFALSE;
  if(mFillHisto) hEvent->Fill(12.5);

  hVtxZ->Fill( vtxPos.z() );

  //ZDC
  StZdcTriggerDetector& ZDC = mMuEvent->zdcTriggerDetector();
  mZDCeast = (UShort_t)ZDC.adcSum(east);
  mZDCwest = (UShort_t)ZDC.adcSum(west);

  //BBC
  StBbcTriggerDetector bbc = mMuEvent->bbcTriggerDetector() ;
  for (UInt_t i = 0; i < bbc.numberOfPMTs(); ++i)
  {
    UInt_t const eastWest = (i < 24) ? 0 : 1 ; // East:0-23, West:24-47
    UInt_t const pmtId    = i % 24 ;         // pmtId:0-23

    if (eastWest == 0) mBbcQ[pmtId] = bbc.adc(i) ;
    else                mBbcQ[pmtId+24] = bbc.adc(i) ;
  }

  Int_t nNodes = mMuDst->numberOfPrimaryTracks();
  if(Debug()){
    LOG_INFO<<"# of primary Tracks in muDst: "<<nNodes<<endm;
  }

  Short_t nTrks    = 0;
  Short_t nBEMCTrks = 0;
  
  //Use default vertex and its tracks.
  //track loop:
  for(Int_t i=0;i<nNodes;i++){
    
    StMuTrack* pMuTrack = mMuDst->primaryTracks(i);
    if(!pMuTrack) continue;

    StMuTrack* gMuTrack = (StMuTrack *)pMuTrack->globalTrack();
    if(!gMuTrack) continue;

    if(!isValidTrack(pMuTrack)) continue;

    mBEMCTraitsIndex[nTrks]  = -999;
    mTPCeTrkFlag[nTrks]      = kFALSE;

    // Calculate global momentum and position at point of DCA to the pVtx
    StThreeVectorF pMom               = pMuTrack->p();
    StPhysicalHelixD gHelix           = gMuTrack->helix(); // Return inner helix (first measured point)
    gHelix.moveOrigin(gHelix.pathLength(vtxPos));
    StThreeVectorF gMom               = gHelix.momentum(mBField*kilogauss);
    StThreeVectorF origin             = gHelix.origin();
  
    mCharge[nTrks]           = pMuTrack->charge();
    mPmag[nTrks]             = pMom.mag();
    mPt[nTrks]               = pMom.perp();
    mEta[nTrks]              = pMom.pseudoRapidity();
    mPhi[nTrks]              = pMom.phi();
    mgPt[nTrks]              = gMom.perp();
    mgEta[nTrks]             = gMom.pseudoRapidity();
    mgPhi[nTrks]             = gMom.phi();
    mgOriginX[nTrks]         = origin.x();
    mgOriginY[nTrks]         = origin.y();
    mgOriginZ[nTrks]         = origin.z();

    mNHitsFit[nTrks]         = pMuTrack->nHitsFit(kTpcId);
    mNHitsPoss[nTrks]        = pMuTrack->nHitsPoss(kTpcId);
    mNHitsDedx[nTrks]        = pMuTrack->nHitsDedx();
    mDedx[nTrks]             = pMuTrack->dEdx()*1.e6; 
    mDndx[nTrks]             = pMuTrack->probPidTraits().dNdxFit();
    mDndxError[nTrks]        = pMuTrack->probPidTraits().dNdxErrorFit();
    mNSigmaE[nTrks]          = pMuTrack->nSigmaElectron();
    mDca[nTrks]              = pMuTrack->dcaGlobal().mag();

    if(mFillHisto){
      hdEdxvsP->Fill(pMom.mag(), mDedx[nTrks]);
      hdNdxvsP->Fill(pMom.mag(), mDndx[nTrks]);
      hnSigEvsP->Fill(pMom.mag(), mNSigmaE[nTrks]);
    }

    //HFT track hits:
    UInt_t  mMap0                     = (UInt_t)(gMuTrack->topologyMap().data(0));
    UChar_t mHftHitsMap               = mMap0>>1 & 0x7F;
    Bool_t  mHasPxl1Hit               = mHftHitsMap>>0 & 0x1;
    Bool_t  mHasPxl2Hit               = mHftHitsMap>>1 & 0x3;
    Bool_t  mHasIstHit                = mHftHitsMap>>3 & 0x3;
    Bool_t  mHasSstHit                = mHftHitsMap>>5 & 0x3;
 
    cout << "Pxl1Hit ~ " << mHasPxl1Hit << endl;
    cout << "Pxl2Hit ~ " << mHasPxl2Hit << endl;
    cout << "IstHit ~ " << mHasIstHit << endl;
    cout << "SstHit ~ " << mHasSstHit << endl;

    //TOF matching:
    mTOFMatchFlag[nTrks] = -1;
    mTOFLocalY[nTrks] = -999.;
    mBeta2TOF[nTrks] = -999.;
    Bool_t matchTofTrack = kFALSE;

    if( &(pMuTrack->btofPidTraits()) ){
      const StMuBTofPidTraits& btofPidTraits = pMuTrack->btofPidTraits();
      mTOFLocalY[nTrks] = btofPidTraits.yLocal();
      mBeta2TOF[nTrks] = btofPidTraits.beta();
      matchTofTrack = btofPidTraits.matchFlag() != 0 ? kTRUE : kFALSE;
      if( matchTofTrack ) mTOFMatchFlag[nTrks] = 1;
      else mTOFMatchFlag[nTrks] = 0;
    }

    //BEMC matching:
    getBemcInfo(pMuTrack,nTrks,nBEMCTrks);
    
    //match BEM:
    if( mBEMCTraitsIndex[nTrks]>=0 || matchTofTrack || mDoMC_ ){nTrks++;}

  }//Track loop

  mNTrks         = nTrks;
  mNBEMCTrks     = nBEMCTrks;

  if(Debug()){
    LOG_INFO<<"# of primary tracks stored: "<<mNTrks<<endm;
    LOG_INFO<<"# of EMC matched Tracks stored: "<<mNBEMCTrks<<endm;
  }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t StUPCTreeMaker::processPicoEvent()
{
  if(mFillHisto) hEvent->Fill(0.5);
  
  StPicoEvent *picoEvent = mPicoDst->event();
  if(!picoEvent){
    LOG_WARN<<"No event level information !"<<endm;
    return kFALSE;
  }
  
  Bool_t validTrigger = kFALSE;

  Int_t nTrigs = 0;
  if(mStreamName.EqualTo("st_upc")){
    for(unsigned i=0;i<mStUPC_TriggerIDs.size();i++){
      if(picoEvent->isTrigger(mStUPC_TriggerIDs[i])){
        validTrigger = kTRUE;
        mTrigId[nTrigs] = mStUPC_TriggerIDs[i];

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
    LOG_WARN<<"No valid UPC related triggers !"<<endm;
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

  hRefMult->Fill( mRefMult );
  
  Int_t nTrks    = 0;

  //track loop:
  for(Int_t i=0;i<nNodes;i++){
    
    StPicoTrack *pTrack = mPicoDst->track(i);
    if(!pTrack) continue;
        
    //track cut
    /*no cut is applied here, cut on analysis level*/
    if(!isValidTrack(pTrack, vtxPos)) continue; 
    
    StThreeVectorF pMom = pTrack->pMom();
    StThreeVectorF gMom = pTrack->gMom();
    StThreeVectorF origin = pTrack->origin();
    
    mCharge[nTrks] = pTrack->charge();
    mPmag[nTrks] = pMom.mag();
    mPt[nTrks] = pMom.perp();
    mEta[nTrks] = pMom.pseudoRapidity();
    mPhi[nTrks] = pMom.phi();
    
    mgPt[nTrks]              = gMom.perp();
    mgEta[nTrks]             = gMom.pseudoRapidity();
    mgPhi[nTrks]             = gMom.phi();
		
    mgOriginX[nTrks]         = origin.x();
		mgOriginY[nTrks]         = origin.y();
		mgOriginZ[nTrks]         = origin.z();

    mNHitsFit[nTrks]         = pTrack->nHitsFit();
    mNHitsPoss[nTrks]        = pTrack->nHitsMax();
    mNHitsDedx[nTrks]        = pTrack->nHitsDedx();
    mDedx[nTrks]             = pTrack->dEdx(); 
    mDndx[nTrks]             = pTrack->dNdx();
    mDndxError[nTrks]        = pTrack->dNdxError();
    mNSigmaE[nTrks]          = pTrack->nSigmaElectron();
    mDca[nTrks]              = (pTrack->dcaPoint()-vtxPos).mag();

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
    mBEMCindex = bemcPidTraitsIndex;

    if( bemcPidTraitsIndex>=0 ){
      StPicoBEmcPidTraits *bemcPidTraits = mPicoDst->bemcPidTraits(bemcPidTraitsIndex);
          
      mBEMCE[nTrks]         = bemcPidTraits->bemcE();
 
    }

    nTrks++;

  }//end of track loop
  

  mNTrks       = nTrks;

  if(Debug()){
    LOG_INFO<<"# of primary tracks stored: "<<mNTrks<<endm;
    //LOG_INFO<<"# of EMC matched Tracks stored: "<<mNBEMCTrks<<endm;
  }
    
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t StUPCTreeMaker::isValidTrack(StMuTrack *pMuTrack) const
{
  Float_t pt  = pMuTrack->pt();
  Float_t eta = pMuTrack->eta(); 
  Float_t dca = pMuTrack->dcaGlobal().mag();

  if(pt<mMinTrkPt)                            return kFALSE;
  if(TMath::Abs(eta)>mMaxTrkEta)              return kFALSE;
  if(pMuTrack->nHitsFit(kTpcId)<mMinNHitsFit) return kFALSE;
  if(pMuTrack->nHitsFit(kTpcId)*1./pMuTrack->nHitsPoss(kTpcId)<mMinNHitsFitRatio)  return kFALSE;
  if(pMuTrack->nHitsDedx()<mMinNHitsDedx)     return kFALSE;
  if(dca>mMaxDca)           return kFALSE;

  return kTRUE;
}
Bool_t StUPCTreeMaker::isValidTrack(StPicoTrack *pTrack, StThreeVectorF vtxPos) const
{
	Float_t pt  = pTrack->pMom().perp();
	Float_t eta = pTrack->pMom().pseudoRapidity();
	Float_t dca = (pTrack->dcaPoint()-vtxPos).mag();

	if(pt<mMinTrkPt)                            return kFALSE;
	if(TMath::Abs(eta)>mMaxTrkEta)              return kFALSE;
	if(pTrack->nHitsFit()<mMinNHitsFit)         return kFALSE;
	if(pTrack->nHitsFit()*1./pTrack->nHitsMax()<mMinNHitsFitRatio) return kFALSE;
	if(pTrack->nHitsDedx()<mMinNHitsDedx)       return kFALSE;
	if(dca>mMaxDca)                             return kFALSE;

	return kTRUE;
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
  mEvtTree->Branch("mBEMCindex", &mBEMCindex, "mBEMCindex/I");
	
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

	//all tracks information
	mEvtTree->Branch("mNTrks", &mNTrks, "mNTrks/I");
  mEvtTree->Branch("mTPCeTrkFlag",mTPCeTrkFlag, "mTPCeTrkFlag[mNTrks]/O");
  mEvtTree->Branch("mBEMCTraitsIndex",mBEMCTraitsIndex,"mBEMCTraitsIndex[mNTrks]/S");

	mEvtTree->Branch("mCharge", mCharge, "mCharge[mNTrks]/I");
  mEvtTree->Branch("mPmag", mPmag, "mPmag[mNTrks]/F");
  mEvtTree->Branch("mPt", mPt, "mPt[mNTrks]/F");
	mEvtTree->Branch("mEta", mEta, "mEta[mNTrks]/F");
  mEvtTree->Branch("mPhi", mPhi, "mPhi[mNTrks]/F");
  mEvtTree->Branch("mgOriginX", mgOriginX, "mgOriginX[mNTrks]/F");
  mEvtTree->Branch("mgOriginY", mgOriginY, "mgOriginY[mNTrks]/F");
  mEvtTree->Branch("mgOriginZ", mgOriginZ, "mgOriginZ[mNTrks]/F");

	mEvtTree->Branch("mgPt", mgPt, "mgPt[mNTrks]/F");
  mEvtTree->Branch("mgEta", mgEta, "mgEta[mNTrks]/F");
  mEvtTree->Branch("mgPhi", mgPhi, "mgPhi[mNTrks]/F");
  
  mEvtTree->Branch("mDedx", mDedx, "mDedx[mNTrks]/F");
  mEvtTree->Branch("mDndx", mDndx, "mDndx[mNTrks]/F");
  mEvtTree->Branch("mNSigmaE", mNSigmaE, "mNSigmaE[mNTrks]/F");
  mEvtTree->Branch("mNHitsFit", mNHitsFit, "mNHitsFit[mNTrks]/I");
	mEvtTree->Branch("mNHitsPoss", mNHitsPoss, "mNHitsPoss[mNTrks]/I");
  mEvtTree->Branch("mNHitsDedx", mNHitsDedx, "mNHitsDedx[mNTrks]/I");
  mEvtTree->Branch("mDndxError", mDndxError, "mDndxError[mNTrks]/F");
	mEvtTree->Branch("mDca", mDca, "mDca[mNTrks]/F");
	
  mEvtTree->Branch("mTOFLocalY", mTOFLocalY, "mTOFLocalY[mNTrks]/F");
	mEvtTree->Branch("mBeta2TOF", mBeta2TOF, "mBeta2TOF[mNTrks]/F");
  mEvtTree->Branch("mTOFMatchFlag", mTOFMatchFlag, "mTOFMatchFlag[mNTrks]/S");

//BEMC pidTrait information
  mEvtTree->Branch("mNBEMCTrks", &mNBEMCTrks, "mNBEMCTrks/S");
  mEvtTree->Branch("mBEMCTrkIndex", mBEMCTrkIndex, "mBEMCTrkIndex[mNBEMCTrks]/S");
  mEvtTree->Branch("mBEMCAdc0", mBEMCAdc0, "mBEMCAdc0[mNBEMCTrks]/S");
  mEvtTree->Branch("mBEMCE0", mBEMCE0, "mBEMCE0[mNBEMCTrks]/F");
  mEvtTree->Branch("mBEMCE", mBEMCE, "mBEMCE[mNBEMCTrks]/F");
  mEvtTree->Branch("mBEMCZDist", mBEMCZDist, "mBEMCZDist[mNBEMCTrks]/F");
  mEvtTree->Branch("mBEMCPhiDist", mBEMCPhiDist, "mBEMCPhiDist[mNBEMCTrks]/F");
  mEvtTree->Branch("mBEMCnEta", mBEMCnEta, "mBEMCnEta[mNBEMCTrks]/F");
  mEvtTree->Branch("mBEMCnPhi", mBEMCnPhi, "mBEMCnPhi[mNBEMCTrks]/F");
}
//_____________________________________________________________________________
void StUPCTreeMaker::bookHistos()
{
	hEvent = new TH1D("hEvent","Event statistics",25,0,25);
  hVtxZ = new TH1D("hVtxZ","hVtxZ",3000,-150,150);
  hNvertex = new TH1D("hNvertex","hNvertex",100,0,100);
  hbestVertex = new TH1D("hbestVertex","hbestVertex",22,-2,20);
  hRefMult = new TH1D("hRefMult","hRefMult",500,0,500);
	hVtxYvsVtxX = new TH2D("hVtxYvsVtxX","hVtxYvsVtxX; V_{x} (cm); V_{y} (cm)",120,-3,3,120,-3,3); 
	hGRefMultvsGRefMultCorr = new TH2D("hGRefMultvsGRefMultCorr","hGRefMultvsGRefMultCorr; grefMultCorr; grefMult",1000,0,1000,1000,0,1000);
    
	hdEdxvsP = new TH2D("hdEdxvsP","hdEdxvsP; p (GeV/c); dE/dx (KeV/cm)",300,0,15,400,0,20);
	hdNdxvsP = new TH2D("hdNdxvsP","hdNdxvsP; p (GeV/c); dN/dx",300,0,15,400,0,200);
	hnSigEvsP = new TH2D("hnSigEvsP","hnSigEvsP; p (GeV/c); n#sigma_{e}",300,0,15,700,-15,20);

}
//_____________________________________________________________________________
void StUPCTreeMaker::printConfig()
{
	const char *decision[2] = {"no","yes"};
	printf("=== Configuration for StUPCTreeMaker ===\n");
	printf("Fill the miniDst tree: %s\n",decision[mFillTree]);
	printf("Fill the QA histo: %s\n",decision[mFillHisto]);
	printf("Maximum |Vz|: %1.2f\n",mMaxVtxZ);
	printf("Minimum Track pt: %1.2f\n",mMinTrkPt);
	printf("Maximum Track |eta| : %1.2f\n",mMaxTrkEta);
	printf("Minimum number of fit hits: %d\n",mMinNHitsFit);
	printf("Minimum ratio of fit hits: %1.2f\n",mMinNHitsFitRatio);
	printf("Minimum number of dedx hits: %d\n",mMinNHitsDedx);
	printf("Maximum dca: %1.2f\n",mMaxDca);
	printf("=======================================\n");
}
//_____________________________________________________________________________
Bool_t StUPCTreeMaker::getBemcInfo(StMuTrack *pMuTrack, const Short_t nTrks, Short_t &nBEMCTrks, Bool_t flag)
{

  Float_t maxtowerE = -999., energy = 0.;
  Float_t zdist = -999., phidist = -999., mindist = 999.;
  Int_t mod = -1, eta=-1, sub=-1;
  Int_t neta = -1, nphi=-1;
  UInt_t maxadc = 0;
  UInt_t maxdsmadc = 0;
  Int_t   softId = -1;

  mEmcCollection = mMuDst->emcCollection();
  if(!mEmcCollection) {
    LOG_WARN << " No Emc Collection for this event " << endm;
    return kFALSE;
  }

  StThreeVectorD position, momentum;
  StThreeVectorD positionBSMDE, momentumBSMDE;
  StThreeVectorD positionBSMDP, momentumBSMDP;

  Double_t bField = mMuDst->event()->runInfo().magneticField()/10.; //Tesla
  Bool_t ok       = kFALSE;
  Bool_t okBSMDE  = kFALSE;
  Bool_t okBSMDP  = kFALSE;
  if(mEmcPosition) {
    ok      = mEmcPosition->projTrack(&position, &momentum, pMuTrack, bField, mEmcGeom[0]->Radius());
    okBSMDE = mEmcPosition->projTrack(&positionBSMDE, &momentumBSMDE, pMuTrack, bField, mEmcGeom[2]->Radius());
    okBSMDP = mEmcPosition->projTrack(&positionBSMDP, &momentumBSMDP, pMuTrack, bField, mEmcGeom[3]->Radius());
  }
  if(!ok) {
    LOG_WARN << " Projection failed for this track ... " << endm;
    return kFALSE;
  }

  Bool_t bemcMatchFlag = kFALSE;
  if(ok && okBSMDE && okBSMDP){

    StSPtrVecEmcPoint& bEmcPoints = mEmcCollection->barrelPoints();
    mindist=1.e9;
    mEmcGeom[0]->getBin(positionBSMDP.phi(), positionBSMDE.pseudoRapidity(), mod, eta, sub); //project on SMD plan
    for(StSPtrVecEmcPointIterator it = bEmcPoints.begin(); it != bEmcPoints.end(); it++) {
      Bool_t associatedPoint = kFALSE;
      StPtrVecEmcCluster& bEmcClusters = (*it)->cluster(kBarrelEmcTowerId);
      if(bEmcClusters.size()==0 ) continue;
      if(bEmcClusters[0]==NULL) continue;
      for(StPtrVecEmcClusterIterator cIter = bEmcClusters.begin(); cIter != bEmcClusters.end(); cIter++){
        Bool_t associatedCluster = kFALSE;
        StPtrVecEmcRawHit& bEmcHits = (*cIter)->hit();
        for(StPtrVecEmcRawHitIterator hIter = bEmcHits.begin(); hIter != bEmcHits.end(); hIter++) {
          if(mod == (Int_t)(*hIter)->module() && eta == (Int_t)(*hIter)->eta() && sub == (Int_t)(*hIter)->sub()) {
            bemcMatchFlag = kTRUE;
            associatedPoint = kTRUE;
            associatedCluster = kTRUE;
            break;
          }
        }
        if(associatedCluster) {
          for(StPtrVecEmcRawHitIterator hitit=bEmcHits.begin(); hitit!=bEmcHits.end();hitit++) {
            if((*hitit)->energy()>maxtowerE) maxtowerE = (*hitit)->energy();
            if((*hitit)->adc()>maxadc) maxadc = (*hitit)->adc();

            softId = (*hitit)->softId(1);
            if(mBemcTriggerSimu && mBemcTriggerSimu->barrelHighTowerAdc(softId)>maxdsmadc) maxdsmadc = mBemcTriggerSimu->barrelHighTowerAdc(softId);
          }
        }
      }

      StPtrVecEmcCluster& smdeClusters = (*it)->cluster(kBarrelSmdEtaStripId);
      StPtrVecEmcCluster& smdpClusters = (*it)->cluster(kBarrelSmdPhiStripId);

      if(associatedPoint) {
        energy += (*it)->energy(); //use point's energy, not tower cluster's energy

        float deltaphi=(*it)->position().phi()-positionBSMDP.phi();
        if(deltaphi>=TMath::Pi()) deltaphi=deltaphi-TMath::TwoPi();
        if(deltaphi<-TMath::Pi()) deltaphi=deltaphi+TMath::TwoPi();

        float rsmdp=mEmcGeom[3]->Radius();
        float pointz=(*it)->position().z();
        float deltaz=pointz-positionBSMDE.z();
        if(sqrt(deltaphi*deltaphi*rsmdp*rsmdp+deltaz*deltaz)<mindist) {
          phidist=deltaphi;
          zdist  =deltaz;
          if(smdeClusters.size()>=1) neta=smdeClusters[0]->nHits();
          if(smdpClusters.size()>=1) nphi=smdpClusters[0]->nHits();
          mindist=sqrt(deltaphi*deltaphi*rsmdp*rsmdp+deltaz*deltaz);
        }
      }//associated
    }
  } // end if (ok && okBSMDE && okBSMDP)

  if(flag) return bemcMatchFlag;

  if(bemcMatchFlag && !flag){
    mBEMCTraitsIndex[nTrks]  = nBEMCTrks;
    mBEMCTrkIndex[nBEMCTrks] = nTrks;
    mBEMCAdc0[nBEMCTrks]     = maxadc;
    mBEMCE0[nBEMCTrks]       = maxtowerE;
    mBEMCE[nBEMCTrks]        = energy;
    mBEMCZDist[nBEMCTrks]    = zdist;  
    mBEMCPhiDist[nBEMCTrks]  = phidist;
    mBEMCnEta[nBEMCTrks]     = neta;
    mBEMCnPhi[nBEMCTrks]     = nphi;
    
    if(Debug()){
      LOG_INFO<<"BEMC associated trkId: "<<pMuTrack->id()<<endm;
      LOG_INFO<<"BEMC associated trkPt: "<<pMuTrack->pt()<<endm;
      LOG_INFO<<"BEMC associated trkEta: "<<pMuTrack->eta()<<endm;
      LOG_INFO<<"BEMC associated trkPhi: "<<pMuTrack->phi()<<endm;
      LOG_INFO<<"BEMC associated trkNHitsFit: "<<pMuTrack->nHitsFit(kTpcId)<<endm;
      LOG_INFO<<"BEMC associated trkNHitsPoss: "<<pMuTrack->nHitsPoss(kTpcId)<<endm;
      LOG_INFO<<"BEMC associated trkNHitsFration: "<<pMuTrack->nHitsFit(kTpcId)*1./pMuTrack->nHitsPoss(kTpcId)<<endm;
      LOG_INFO<<"BEMC associated trkNHitsDedx: "<<pMuTrack->nHitsDedx()<<endm;
      LOG_INFO<<"BEMC associated trkDca: "<<pMuTrack->dcaGlobal().mag()<<endm;
      LOG_INFO<<"BEMC Adc0: "<<mBEMCAdc0[nBEMCTrks]<<endm;
      LOG_INFO<<"BEMC E0: "<<mBEMCE0[nBEMCTrks]<<endm;
      LOG_INFO<<"BEMC E: "<<mBEMCE[nBEMCTrks]<<endm;
      LOG_INFO<<"BEMC ZDist: "<<mBEMCZDist[nBEMCTrks]<<endm;
      LOG_INFO<<"BEMC PhiDist: "<<mBEMCPhiDist[nBEMCTrks]<<endm;
      LOG_INFO<<"BEMC nEta: "<<(Int_t)mBEMCnEta[nBEMCTrks]<<endm;
      LOG_INFO<<"BEMC nPhi: "<<(Int_t)mBEMCnPhi[nBEMCTrks]<<endm;
    }

    nBEMCTrks++;
  };

  return bemcMatchFlag;
}
//____________________________________________________________________________

