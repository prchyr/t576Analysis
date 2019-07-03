/*
  copyright s. prohira and the T576 collaboration 2018
  released under the GNU General Public License version 3
*/


/*
This is a storage class for useful properties of events. it then allows for access via ttree->Draw(), etc. 

 */

#ifndef TEVENT_SUMMARY_H
#define TEVENT_SUMMARY_H


#include "TROOT.h"
#include "TRint.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TObject.h"
#include "TNamed.h"
#include "TLine.h"
#include "TFile.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "cnpy.h"
#include "TSystemFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "TVirtualFFT.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"
#include "TVectorD.h"
#include "TMatrixTBase.h"
#include "TMatrixD.h"
#include "TVector3.h"
#include "TDecompSVD.h"
#include <iostream>
#include <fstream>
#include <vector>




using namespace std;





class TEventSummary : public TObject{
public:
  /************constructors**********/
  TEventSummary(int dataset=0){
    fDATASET=dataset;
  }
  ~TEventSummary(){
    //delete []fSurfData;
  }

  /*
    which dataset:
    0-raw data
    1-SVD filtered
    2-null data
    3-null SVD filtered

  */
  int dataset;
  TVector3 txPos;
  double txAng, txDist;
  ULong64_t timestamp, scopeTime, surfTime;
  double frequency;
  double power;
  double analogBandwidth;
  /*the type of antenna for:
    antennaType[0] = TX, 
    antennaType[1] = scope ch1
    antennaType[2] = scope ch2
    antennaType[3] = scope ch3
    antennaType[4] = scope ch4 (always ICT see below)
    antennaType[5] = surf (always =lpda, see below)

    the key is:
    1=vivaldi
    2=wifi
    3=lpda
    4=Jiwoo (lpda with a dish)
    5=ICT
  */
  double antennaType[6]={0.,0.,0.,0.,0.,0.};
  //various event numbering indices.
  int subEvNo, scopeEvNo, surfEvNo, evNo, surfNEvents, scopeNEvents;
  //accesors for the file names associated with the event.
  TString * scopeFilename=new TString();
  TString * surfFilename=new TString();

  //static coordinates useful for stuff

  TVector3 beamPipeExit;//(.6, 0., -3.6);
  TVector3 targetFront;//(.6, 0., -2.);
  TVector3 targetBack;//(.6, 0., 2.);
  TVector3 beamDump;//(.6, 0., 4.);
  TVector3 backWall;//(.6, 0., 5.6);
  
  //useful flags and info
  //run major and minor
  int major=0, minor=0;
  //major and minor for the background run associated with this run
  int backmajor=0, backminor=0;
  //is the tx on
  int txOn=0;
  //is the surf on
  int surfOn=0;
  //polarization (0=V, 1=H)
  int polarization=0;
  //is the event good?
  int isGood=0;

  double charge=0;

  /* some useful power and field quantities */

  double pow_0_10[3]={0.,0.,0.};
  double pow_10_20[3]={0.,0.,0.};
  double pow_20_30[3]={0.,0.,0.};
  double pow_30_40[3]={0.,0.,0.};
  double pow_40_50[3]={0.,0.,0.};
  double pow_50_60[3]={0.,0.,0.};
  double pow_60_70[3]={0.,0.,0.};
  double pow_70_80[3]={0.,0.,0.};
  double pow_80_90[3]={0.,0.,0.};
  double pow_90_100[3]={0.,0.,0.};

  double rms_0_10[3]={0.,0.,0.};
  double rms_10_20[3]={0.,0.,0.};
  double rms_20_30[3]={0.,0.,0.};
  double rms_30_40[3]={0.,0.,0.};
  double rms_40_50[3]={0.,0.,0.};
  double rms_50_60[3]={0.,0.,0.};
  double rms_60_70[3]={0.,0.,0.};
  double rms_70_80[3]={0.,0.,0.};
  double rms_80_90[3]={0.,0.,0.};
  double rms_90_100[3]={0.,0.,0.};

  double peakPower[3]={0.,0.,0.};
  double peakV[3]={0.,0.,0.};
  double tOfPeakV[3]={0.,0.,0.};
  double peakHilbert[3]={0.,0.,0.};
  double tOfPeakHilbert[3]={0.,0.,0.};

  TVector3 pos[3];
  //unused currently
  double cableLengths[3]={0.,0.,0.};
  double velocityFactor=.82;
  
private:
  int fNEntriesSurf=0, fNEntriesScope=0;
  double fInterpGSs=0.;
  TString fScopeFilename="20181025191401run0_4.root", fInstallDir="20181025191401run0_4.root";
  TString fSurfFilename="20181025191401run0_4.py";
  int fDATASET=0;
  ClassDefNV(TEventSummary, 1)
};

#endif

  

  

