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
//#include "TIter.h"
#include "TSystemFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"
#include <iostream>
#include <fstream>
#include <vector>

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"


using namespace CLHEP;
using namespace std;


class T576Event : public TObject{
public:
  /************constructors**********/
  T576Event(){checkStatus();}
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
  Hep3Vector txPos;
  double txAng, txDist;
  double charge;
  ULong64_t timestamp;
  double frequency;
  double power;
  double antennaType[6]={0.,0.,0.,0.,0.,0.};
  //various event numbering indices.
  int subEvNo, scopeEvNo, surfEvNo, evNo;
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
  int loadScopeEvent(int run_major, int run_minor, int event);
  //load an event using a running scope event number index. this index is built the first time the T576Event class is used, starting at 0 for the first event in run 0_0 and going upward.
  //when you use this, the correct positions and stuff will be loaded for the correct run etc.
  int loadScopeEvent(int event);
  //build the event index ROOT file which allows for rapid grabbing of the right data (done once, first time you use the class. you can also force a re-build if you change the data directory [don't.]) 
  int buildEventIndex(int force=0);
  //run on construction to make sure the indexing is valid and useable.
  int checkStatus();
  //set the interpolation level
  void setInterpGsNs(double value){fInterpGsNs=value;};

  /**************utility and analysis functions**************/

  //get the charge in this event from the ict trace.
  double getCharge(TGraph * ict);
  //provide a graph and a pointer to a graph, and the desired samplerate in GS/s
  int getInterpolatedGraph(TGraph * inGraph, TGraph *outGraph);
  
private:

  double fInterpGsNs=0.;
  TString fScopeFilename, fInstallDir;
  int fIndexBuilt=0;
  TTreeIndex * fScopeEvNoIndex, *fMajorMinorIndex, *fRunLogIndex;
  TTree *fIndexTree=new TTree();
  TTree *fEventTree=new TTree();
  TFile *fEventFile=new TFile();
  TFile * fRunLog= new TFile();
  TTree * fRunLogTree = new TTree();
  
public:
  class Scope {
  public:
    //receiver positions
    Hep3Vector pos[4];
    double dist[4], ang[4];
    //double arrays for the individual traces
    double  ch[4][20000];
    //and the x axis time of the recorded traces
    double  time[20000];
    //tgraphs of each event channel
    TGraph * gr[4]={new TGraph(), new TGraph(),new TGraph(),new TGraph()};


    //loads the antenna position iformation into the pos[] above, using the run log.
    int getAntennaPositions(int run_major, int run_minor);
  private:
    
    ClassDefNV(Scope, 1);
  };


  class Surf{
  public:
    Hep3Vector pos[12];
    double dist[12], ang[12];
    double  ch[12][1024];
    TGraph * gr[12]={new TGraph()};
    double  time[1024];

    TGraph2D * map=0;
    int getAntennaPositions(int run_major, int run_minor);
    TGraph2D *buildMap(int mmStep=500);

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
