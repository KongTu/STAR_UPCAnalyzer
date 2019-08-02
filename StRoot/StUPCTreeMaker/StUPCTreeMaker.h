#ifndef STUPCTREEMAKER_HH
#define STUPCTREEMAKER_HH

/***************************************************************************
 *
 * $Id: StUPCTreeMaker.h 2015/04/09  Exp $ 
 * StUPCTreeMaker - class to produce miniTree for UPC related analysis
 * Author: Kong
 *--------------------------------------------------------------------------
 *
 ***************************************************************************/
#include "StMaker.h"
#include "StTreeStructure.h"

#include "StThreeVectorF.hh"
#include "TLorentzVector.h"
#include "TComplex.h"

#include "StMtdUtil/StMtdConstants.h"

#include <vector>
#include <map>
#ifndef ST_NO_NAMESPACES
using std::vector;
#endif

class TH1D;
class TH2D;
class TString;
class TTree;
class TFile;

class StMuDstMaker;
class StMuDst;
class StEmcCollection;
class StEmcPosition;
class StEmcGeom;
class StEmcRawHit;
class StMuTrack;
class StMuMtdHit;

class StPicoDstMaker;
class StPicoDst;
class StPicoTrack;
class StPicoMtdHit;

#if !defined(ST_NO_TEMPLATE_DEF_ARGS) || defined(__CINT__)
typedef vector<Int_t> IntVec;
typedef vector<Double_t> DoubleVec;
typedef vector<TLorentzVector> LorentzVec;
#else
typedef vector<Int_t, allocator<Int_t>> IntVec;
typedef vector<Double_t, allocator<Double_t>> DoubleVec;
typedef vector<TLorentzVector, allocator<TLorentzVector>> LorentzVec;
#endif

const Int_t          NPCAbin=10;

class StUPCTreeMaker : public StMaker {
	public:
		StUPCTreeMaker(const Char_t *name = "StUPCTreeMaker");
		~StUPCTreeMaker();

		Int_t    Init();
		Int_t    InitRun(const Int_t runNumber);
		Int_t    Make();
		Int_t    Finish();

		void     setMaxVtxR(const Double_t max);
		void     setMaxVtxZ(const Double_t max);
		void     setMaxVzDiff(const Double_t max);
		void     setMinTrackPt(const Double_t min);
		void     setMaxTrackEta(const Double_t max);
		void     setMinNHitsFit(const Int_t min);
		void     setMinNHitsFitRatio(const Double_t min);
		void     setMinNHitsDedx(const Int_t min);
		void     setMaxDca(const Double_t max);
		void     setMaxnSigmaE(const Double_t max);
		void     setMaxBeta2TOF(const Double_t max);
		void     setFillHisto(const Bool_t fill);
		void     setFillTree(const Bool_t fill);
		void     setOutFileName(const TString name);
		void     setStreamName(const TString name);
		void     setPrintMemory(const Bool_t pMem);
		void     setPrintCpu(const Bool_t pCpu);
		void     setPrintConfig(const Bool_t print);
		void	 setMCevent(const Bool_t isMC);

	protected:
		void     printConfig();
		void     bookTree();
		void     bookHistos();
		Bool_t   processMuDstEvent();
		Bool_t   processPicoEvent();
		void     fillEventPlane();
		Bool_t   isValidTrack(StMuTrack *pMuTrack) const;
		Bool_t   isValidTrack(StPicoTrack *pTrack, StThreeVectorF vtxPos) const;

		void     initEmc();
		void     finishEmc();
		void     buildEmcIndex();

		bool     getBEMC(const StMuTrack* t, int* id, int* adc, float* ene, float* d, int* nep, int* towid);
		Bool_t   getBemcInfo(StMuTrack *pMuTrack, const Short_t nTrks, Short_t &nEMCTrks, Bool_t flag=kFALSE);
		Float_t  getBemcInfo(StMuTrack *pMuTrack);
		Bool_t   isMtdHitFiredTrigger(const StPicoMtdHit *hit);
		Bool_t   isMtdHitFiredTrigger(const StMuMtdHit *hit);
		Bool_t   isQTFiredTrigger(const Int_t qt, const Int_t pos);
		void     triggerData();


	private:

		StMuDstMaker    *mMuDstMaker;          // Pointer to StMuDstMaker
		StMuDst         *mMuDst;              // Pointer to MuDst event
	
		StPicoDstMaker  *mPicoDstMaker;
		StPicoDst       *mPicoDst;

		StTriggerSimuMaker* mTriggerSimuMaker;
        StBemcTriggerSimu*  mBemcTriggerSimu;

		//variable for tree
		Int_t    mRunId;
		Int_t    mEventId;
		Int_t    mNTrigs;
		Int_t    mTrigId[64];
		Int_t    mBEMCindex;
		Short_t  mRefMult;
		Short_t  mGRefMult;

		Float_t  mBBCRate;
		Float_t  mZDCRate;
		Float_t  mBField;
		Float_t  mVpdVz;
		Float_t  mVertexX;
		Float_t  mVertexY;
		Float_t  mVertexZ;

		//BBC
		Int_t    mBbcQ[48];

		//nTOF
		Int_t    mNTofHits;

		//ZDC
		Int_t mZDCwest;
		Int_t mZDCeast;

		//track information
		Int_t  mNTrks;
		Int_t   mCharge[mMax];
		Float_t  mPmag[mMax];
		Float_t  mPt[mMax];
		Float_t  mEta[mMax];
		Float_t  mPhi[mMax];

		Float_t  mgPt[mMax];
		Float_t  mgEta[mMax];
		Float_t  mgPhi[mMax];

		Float_t  mgOriginX[mMax];
		Float_t  mgOriginY[mMax];
		Float_t  mgOriginZ[mMax];
		
		Int_t mNHitsFit[mMax];
		Int_t mNHitsPoss[mMax];
		Int_t mNHitsDedx[mMax];

		Float_t mDedx[mMax];
		Float_t mDndx[mMax];
		Float_t mDndxError[mMax];
		Float_t mNSigmaE[mMax];
		Float_t mDca[mMax];

		Float_t  mTOFLocalY[mMax];
		Short_t  mTOFMatchFlag[mMax];
		Float_t  mBeta2TOF[mMax];
		Bool_t   mTPCeTrkFlag[mMax];
		Short_t  mBEMCTraitsIndex[mMax];
		
		//BEMC pidTrait information
		Short_t  mNBEMCTrks;
		Short_t  mBEMCTrkIndex[mMax];
		Short_t  mBEMCAdc0[mMax];
		Float_t  mBEMCE0[mMax];
		Float_t  mBEMCE[mMax];
		Float_t  mBEMCZDist[mMax];
		Float_t  mBEMCPhiDist[mMax];
		Float_t  mBEMCnEta[mMax];
		Float_t  mBEMCnPhi[mMax];

		//StRefMultCorr *refMultCorr; //decide centrality
		Bool_t         mFillTree;            // Flag of fill the event tree
		Bool_t         mFillHisto;           // Flag of fill the histogram
		Bool_t         mPrintConfig;         // Flag to print out task configuration
		Bool_t         mPrintMemory;         // Flag to print out memory usage
		Bool_t         mPrintCpu;            // Flag to print out CPU usage
		Bool_t		   mDoMC_;
		TString        mStreamName;          // Data stream name
		TFile          *fOutFile;            // Output file
		TString        mOutFileName;         // Name of the output file 
		StEvtData      mEvtData;
		TTree          *mEvtTree;            // Pointer to the event tree

		Double_t       mMaxVtxR;             // Maximum vertex r
		Double_t       mMaxVtxZ;             // Maximum vertex z
		Double_t       mMaxVzDiff;           // Maximum VpdVz-TpcVz
		Double_t       mMinTrkPt;            // Minimum track pt
		Double_t       mMaxTrkEta;           // Maximum track eta
		Int_t          mMinNHitsFit;         // Minimum number of hits used for track fit
		Double_t       mMinNHitsFitRatio;    // Minimum ratio of hits used for track fit
		Int_t          mMinNHitsDedx;        // Minimum number of hits used for de/dx
		Double_t       mMaxDca;              // Maximum track dca
		Double_t       mMaxnSigmaE;          // Maximum nSigmaE cut
		Double_t       mMaxBeta2TOF;         // Maximum |1-1./beta| for TpcE
		
		StEmcCollection *mEmcCollection;
		StEmcPosition   *mEmcPosition;
		StEmcGeom       *mEmcGeom[4];
		StEmcRawHit     *mEmcIndex[4800];

		IntVec         mStPhysics_TriggerIDs;
		IntVec 		   mStUPC_TriggerIDs;

		//define histograms ongoing...
		TH1D           *hEvent;
		TH1D		   *hVtxZ;
		TH1D		   *hNvertex;
		TH1D		   *hbestVertex;
		TH1D		   *hRefMult;
		TH2D           *hVtxYvsVtxX;
		TH2D           *hVPDVzvsTPCVz;
		TH1D           *hVzDiff;
		TH2D           *hGRefMultvsGRefMultCorr;
		TH1D           *hCentrality;

		TH2D           *hdEdxvsP;
		TH2D           *hdNdxvsP;
		TH2D           *hnSigEvsP;
		TH2D           *hBetavsP;


		ClassDef(StUPCTreeMaker, 1)
};
inline void StUPCTreeMaker::setMaxVtxR(const Double_t max) { mMaxVtxR = max; }
inline void StUPCTreeMaker::setMaxVtxZ(const Double_t max) { mMaxVtxZ = max; }
inline void StUPCTreeMaker::setMaxVzDiff(const Double_t max) { mMaxVzDiff = max; }
inline void StUPCTreeMaker::setMinTrackPt(const Double_t min){ mMinTrkPt = min;}
inline void StUPCTreeMaker::setMaxTrackEta(const Double_t max){ mMaxTrkEta = max; }
inline void StUPCTreeMaker::setMinNHitsFit(const Int_t min) { mMinNHitsFit = min; }
inline void StUPCTreeMaker::setMinNHitsFitRatio(const Double_t min) { mMinNHitsFitRatio = min; }
inline void StUPCTreeMaker::setMinNHitsDedx(const Int_t min) { mMinNHitsDedx = min; }
inline void StUPCTreeMaker::setMaxDca(const Double_t max) { mMaxDca = max; }
inline void StUPCTreeMaker::setMaxnSigmaE(const Double_t max) { mMaxnSigmaE = max; }
inline void StUPCTreeMaker::setMaxBeta2TOF(const Double_t max) { mMaxBeta2TOF = max; }
inline void StUPCTreeMaker::setFillHisto(const Bool_t fill) { mFillHisto = fill; }
inline void StUPCTreeMaker::setFillTree(const Bool_t fill) { mFillTree = fill; }
inline void StUPCTreeMaker::setOutFileName(const TString name) { mOutFileName = name; }
inline void StUPCTreeMaker::setStreamName(const TString name) { mStreamName = name; }
inline void StUPCTreeMaker::setPrintMemory(const Bool_t pMem) { mPrintMemory = pMem; }
inline void StUPCTreeMaker::setPrintCpu(const Bool_t pCpu) { mPrintCpu = pCpu; }
inline void StUPCTreeMaker::setPrintConfig(const Bool_t print) { mPrintConfig = print; }
inline void StUPCTreeMaker::setMCevent(const Bool_t isMC){ mDoMC_ = isMC; }
#endif
