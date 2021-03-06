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
#include "TNtuple.h"
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
  T576Event(double interpGSs=0., int useFilteredData=0){
    checkStatus();
    fUSE_FILTERED_DATA=useFilteredData;
    if (interpGSs!=0){
      setInterpGSs(interpGSs);
    }
  }
  //load a T576 event using run major and minor. actual file name not needed.
  T576Event(int run_major, int run_minor, int event,int useFilteredData=0){
    major=run_major;
    minor=run_minor;
    checkStatus();
    fUSE_FILTERED_DATA=useFilteredData;
    loadScopeEvent(run_major, run_minor,event);
  }
  ~T576Event(){
    //delete []fSurfData;
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
  TVector3 targetCenter;//(.6, 0., 0.);
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

  /****************operational functions*****************/
  
  //load an event from a particular run file. actual file name not needed, just the major and minor associated with it. so XXXXrun4_6.root event number 15 is loaded by doing loadScopeEvent(4, 6, 15); 
  int loadScopeEvent(int run_major, int run_minor, int event, bool remove_dc_offset=true);
  //same for surf
  int loadSurfEvent(int run_major, int run_minor, int event, bool remove_dc_offset=true);
  //load an event using a running scope event number index. this index is built the first time the T576Event class is used, starting at 0 for the first event in run 0_0 and going upward.
  //when you use this, the correct positions and stuff will be loaded for the correct run etc.
  int loadScopeEvent(int event, bool remove_dc_offset=true);
  //same for surf
  int loadSurfEvent(int event, bool remove_dc_offset=true);
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
  //make a pointing map from the SURF. several ways:
  //type==0 is using the raw voltage versus time trace
  //type==1 uses the hilbert envelope of the voltage versus time trace
  //type==2 uses the power versus time trace
  TH2D* pointingMap(double dx=.3, int draw=1, int type=1);
  TH2D* pointingMapDev(double dx, int draw, int hilbert, TVector3 *position);
  TVector3 * fixPositionsDev(double dx, int maxIter, int hilbert, TVector3 source, TVector3 *positions);
  vector<TPolyLine*> getTheRoom(Color_t lineColor=kRed, Color_t tgtColor=kGreen);
  TGraph * getAntennas(Color_t color=kBlue);
  double getInterpGSs();

    //draw a bunch of graphs from a t576event object
  //the align parameter should be 0 for non-alignement, otherwise it is
  //the number of ns of allowed wiggle.
  vector<TGraph *> draw(int major, int minor, int scopeOrSurf, int channel, int num, double align=5., double tLow=0., double tHigh=999999., TString drawOption="al PLC");
   //draw average of a bunch of graphs from a t576event object
  TGraph * drawAvg(int major, int minor, int scopeOrSurf, int channel, int num, double align=5., double tLow=0., double tHigh=999999.,TString drawOption="al PLC");
  TGraph * drawAll(int major, int minor, int scopeOrSurf, int channel, int num, double align=5., double tLow=0., double tHigh=999999., TString drawOption="al PLC");
  TGraph * drawAvgHilbert(int major, int minor, int scopeOrSurf, int channel, int num, double align=5., double tLow=0., double tHigh=999999.,TString drawOption="al PLC");
    TGraph * avgHilbert(int major, int minor, int scopeOrSurf, int channel, int num, double align=5., double tLow=0., double tHigh=999999.);
  //draw an average of a bunch of spectrograms of t576event object
  TH2D * drawAvgSpectrogram(int major, int minor, int scopeOrSurf, int channel, int num, Int_t binsize , Int_t overlap, Int_t zero_pad_length, int win_type, int dbFlag, int start=0);
  TH2D * drawAvgSpectrogram(int major, int minor, int scopeOrSurf,int channel, int num, double tLow, double tHigh, Int_t binsize , Int_t overlap, Int_t zero_pad_length, int win_type, int dbFlag, int start=0);
  TGraph * drawAvgPSD(int major, int minor, int scopeOrSurf,int channel, int num, double tLow, double tHigh, int dbFlag);
  TGraph * drawAvgPSD(int major, int minor, int scopeOrSurf,int channel, int num, int dbFlag);

  TNtuple * integrateAllWithSideband(int major, int minor, int scopeOrSurf, int channel, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, double sbxmin, double sbxmax, double sbymin, double sbymax, double scale=1.);
  
  TNtuple * sidebandSubtractAll(int major, int minor, int scopeOrSurf, int channel, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, double scale=1.);


  
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
  cnpy::NpyArray fTimesArray;
  double * fSurfTimes;//=(double*) malloc(20000);
  TFile * fRunLog= new TFile();
  TTree * fRunLogTree = new TTree();

  bool fScopeLoaded=false;
  bool fSurfLoaded=false;
  int fUSE_FILTERED_DATA=0;
  
public:
  class Scope  {
  public:
    //    Scope(): fTxPos(txPos){};
    //receiver positions
    double cableLengths[4]={0., 0., 0., 10.};
    double velocityFactor=-.82;
    double delays[4];
    
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
