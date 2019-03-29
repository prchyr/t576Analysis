/*
copyright s. prohira and the T576 collaboration 2018
released under the GNU General Public License version 3
*/

#ifndef T576EVENT_BASE_H
#define T576EVENT_BASE_H

#include "TROOT.h"
#include "TRint.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TObject.h"
#include "TFile.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "cnpy.h"
//#include "TIter.h"
#include "TSystemFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TVector3.h"
#include "TMath.h"
#include "TColor.h"
#include "TStyle.h"
#include "TPolyLine.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <vector>




 #include "TUtil.hh"


using namespace std;
using namespace TUtil;

class T576Event : public TObject{
public:
  /************constructors**********/
  T576Event(double interpGSs=0.){
    checkStatus();
    if (interpGSs!=0){
      setInterpGSs(interpGSs);
    }
  }
  //load a T576 event using run major and minor. actual file name not needed.
  T576Event(int run_major, int run_minor, int event){
    major=run_major;
    minor=run_minor;
    checkStatus();
    loadScopeEvent(run_major, run_minor,event);
  }

  /*************globals for every event*****************/
  //the two DAQs. 
  class Scope;
  class Surf;

  //useful variables, gathered from the run log for each event.
  TVector3 txPos;
  double txAng, txDist;
  double charge;
  ULong64_t timestamp, scopeTime, surfTime;
  double frequency;
  double power;
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

  //useful flags and info
  //run major and minor
  int major=0, minor=0;
  //is the tx on
  int txOn=0;
  //is the surf on
  int surfOn=0;
  //polarization (0=V, 1=H)
  int polarization=0;
  //is the event good?
  int isGood=0;

  /****************operational functions*****************/
  
  //load an event from a particular run file. actual file name not needed, just the major and minor associated with it. so XXXXrun4_6.root event number 15 is loaded by doing loadScopeEvent(4, 6, 15); 
  int loadScopeEvent(int run_major, int run_minor, int event, bool remove_dc_offset=true);
  //same for surf
  int loadSurfEvent(int run_major, int run_minor, int event);
  //load an event using a running scope event number index. this index is built the first time the T576Event class is used, starting at 0 for the first event in run 0_0 and going upward.
  //when you use this, the correct positions and stuff will be loaded for the correct run etc.
  int loadScopeEvent(int event, bool remove_dc_offset=true);
  //same for surf
  int loadSurfEvent(int event);
  //build the event index ROOT file which allows for rapid grabbing of the right data (done once, first time you use the class. you can also force a re-build if you change the data directory [don't.]) 
  int buildEventIndex(int force=0);
  //run on construction to make sure the indexing is valid and useable.
  TString getSurfFilename(int inMaj, int inMin);
  //called on construction, checks that things are where they should be.
  int checkStatus();
  //set the interpolation level. this is used in the graphs, but the
  //channel arrays are always the raw values.
  void setInterpGSs(double value){fInterpGSs=value;};

  /**************utility and analysis functions**************/

  //get the charge in this event from the ict trace.
  double getCharge(TGraph * ict);
  int drawGeom(int scopeChan=999, int surfChan=999);
  TH2D* pointingMap(double dx=.3, int draw=1, int hilbert=1);
  TH2D* pointingMapDev(double dx, int draw, int hilbert, TVector3 *position);
  TVector3 * fixPositionsDev(double dx, int maxIter, int hilbert, TVector3 source, TVector3 *positions);
  vector<TPolyLine*> getTheRoom();
  TGraph * getAntennas();

  
private:
  int fNEntriesSurf=0, fNEntriesScope=0;
  double fInterpGSs=0.;
  TString fScopeFilename="20181025191401run0_4.root", fInstallDir="20181025191401run0_4.root";
  TString fSurfFilename="20181025191401run0_4.py";
  int fIndexBuilt=0;
  TTreeIndex * fScopeEvNoIndex=new TTreeIndex();
  TTreeIndex * fSurfEvNoIndex=new TTreeIndex();
  TTreeIndex *fMajorMinorIndex=new TTreeIndex();
  TTreeIndex * fScopeTimeIndex=new TTreeIndex();
  TTreeIndex * fSurfTimeIndex=new TTreeIndex();
  TTreeIndex *fRunLogIndex=new TTreeIndex();
  TTree *fIndexTree=new TTree();
  TTree *fEventTree=new TTree();
  TFile *fEventFile=new TFile();
  cnpy::npz_t fDataset;
  cnpy::NpyArray fSurfDataArray;
  short * fSurfData;//=(short*) calloc(60000000, sizeof(short));
  double * fSurfTimes;//=(double*) malloc(20000);
  TFile * fRunLog= new TFile();
  TTree * fRunLogTree = new TTree();

  bool fScopeLoaded=false;
  bool fSurfLoaded=false;
  
public:
  class Scope  {
  public:
    //    Scope(): fTxPos(txPos){};
    //receiver positions
    TVector3 pos[4];
    double dist[4], ang[4];
    //double arrays for the individual traces
    double  dat[4][20000];
    //and the x axis time of the recorded traces
    double  time[20000];
    //tgraphs of each event channel
    TGraph * ch[4]={new TGraph(), new TGraph(),new TGraph(),new TGraph()};
    void draw();
    

  private:
    
    ClassDefNV(Scope, 1);
  };


  class Surf{
  public:
    Surf(){
      double incr=1./3.2;
      for(int i=0;i<1024;i++){
	time[i]=i*incr;
      }
    }
    //these are measured in feet.
    double cableLengths[12]={50., 20., 20., 10., 0., 0., 0., 0., 0., 0., 10., 20.};
    double velocityFactor=-.82;//from Krijn and Kian, negative for math reasons
    
    TVector3 pos[12];
    double dist[12], ang[12];
    double delays[12];
    double  dat[12][1024];
    TGraph * ch[12]={new TGraph(), new TGraph(),new TGraph(),new TGraph(),new TGraph(),new TGraph(),new TGraph(),new TGraph(),new TGraph(),new TGraph(),new TGraph(),new TGraph()};
    double  time[1024];

    TGraph2D * map=0;
    int drawGeom();
    void draw();
    //int drawMap(int mStep=.5);

  private:

    int fMapMade=0;

    
    ClassDefNV(Surf, 1);
  };

  //public things
  Surf * surf=new Surf();
  Scope * scope=new Scope();


  
  ClassDefNV(T576Event, 1);
};

#endif
