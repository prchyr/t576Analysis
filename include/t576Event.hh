#ifndef T576EVENT_BASE_H
#define T576EVENT_BASE_H

#include "TROOT.h"
#include "TRint.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraph2D.h"

#include <vector>

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"


using namespace CLHEP;
using namespace std;


class t576Event : public TObject{
public:
  //constructors
  t576Event(){}
  //load a t576 event
  t576Event(int run_major, int run_minor, int event):loadEvent(run_major, run_minor){}
  //load a t576 event from a radio scatter event
  t576Event(RadioScatterEvent *rs):loadEvent(rs){}
  //load a t576 event from the generic tree in the raw capture files
  t576Event(TTree * tree):loadEvent(tree){}


  //global info for all channels in an event.
  double charge;
  long timestamp;
  double frequency;
  double power;
  Hep3Vector txPos;
  int major, minor;


  //the two DAQs. 
  class Scope;
  class Surf;
private:
  int fRunLoaded=0;

  class Scope {
  public:
    Hep3Vector pos[4]={(0.,0.,0.)};
    TGraph * ch[4]={0};
    double * time;
    ClassDefNV(Scope, 1);
  };


  class Surf{
  public:
    Hep3Vector pos[12]={(0.,0.,0.)};
    TGraph * ch[12]={0};
    double * time;

    TGraph2D * map=0;
  private:
    TGraph2D *buildMap(int mmStep=500);
    int fMapMade=0;

    
    ClassDefNV(Surf, 1);
  };
  ClassDefNV(t576Event, 1);
}
