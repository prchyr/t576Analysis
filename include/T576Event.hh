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
  //constructors
  T576Event(){checkStatus();}
  //load a T576 event
  //  T576Event(int run_major, int run_minor, int event){;}
  T576Event(int run_major, int run_minor, int event){
    major=run_major;
    minor=run_minor;
    checkStatus();
    loadEvent(run_major, run_minor,event);
  }
  //load a T576 event from a radio scatter event
  //T576Event(int run_major, int run_minor,RadioScatterEvent *rs):loadEvent(run_major, run_minor, rs){}
  //load a T576 event from the generic tree in the raw capture files
  T576Event(int run_major, int run_minor,TTree * inTree){;}//{loadEvent(run_major, run_minor, tree);}



  //the two DAQs. 
  class Scope;
  class Surf;

  //global info for all channels in an event.
  double charge;
  ULong64_t timestamp;
  double frequency;
  double power;
  int subEvNo, scopeEvNo, surfEvNo, evNo;
  TString * scopeFilename=new TString();
  TString * surfFilename=new TString();
  Hep3Vector txPos;
  int major, minor, event;

  int loadEvent(int run_major, int run_minor, int event);
  //  void loadEvent(int run_major, int run_minor, RadioScatterEvent *rs);
  int loadEvent(int run_major, int run_minor, TTree *inTree);
  int getCharge(TGraph * ict);
  int buildEventIndex(int force=0);
  int checkStatus();
  
private:
  int fRunLoaded=0;
  TString fScopeFilename, fInstallDir;
  int fIndexBuilt=0;
  TTreeIndex * fScopeEvNoIndex, *fMajorMinorIndex;
  TTree *fIndexTree=new TTree();
  
public:
  class Scope {
  public:
    Hep3Vector pos[4];
    double  ch[4][20000];
    //    TGraph *fGr=new TGraph();
    TGraph * gr[4]={new TGraph(), new TGraph(),new TGraph(),new TGraph()};
    double  time[20000];
    int getAntennaPositions(int run_major, int run_minor);
  private:

    ClassDefNV(Scope, 1);
  };


  class Surf{
  public:
    Hep3Vector pos[12];
    double  ch[12][1024];
    TGraph * gr[12]={new TGraph()};
    double  time[1024];

    TGraph2D * map=0;
    int getAntennaPositions(int run_major, int run_minor);
  private:
    TGraph2D *buildMap(int mmStep=500);
    int fMapMade=0;

    
    ClassDefNV(Surf, 1);
  };

  //public things
  Surf surf;
  Scope scope;


  
  ClassDefNV(T576Event, 1);
};

#endif
