/*
  copyright s. prohira and the T576 collaboration 2018
  released under the GNU General Public License version 3
*/

/*
This is a general analysis class for the T576 data. It relies on a dataset comprised of T576Events. The primitive methods are pulled from TUtil.hh.

 */

#ifndef TANALYZE_BASE_H
#define TANALYZE_BASE_H


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



#include "TUtilGraph.hh"
#include "TUtil.hh"
#include "T576Event.hh"


using namespace std;


namespace TAnalyze{
  //some drawing functions

  TH2D* avgSpectrogram(int ch, int dataset, int major, int minor, int nfft=512, int nOverlap=450, int padTo=800, int window=3, int log=0, int norm=0);
  TH2D* avgSpectrogram(int ch, int dataset, int major, int minor, double tLow, double tHigh, int nfft, int nOverlap, int padTo, int window, int log, int norm);
  //average spectrogram from a run compared to null data
  int drawAvgRealNull(int ch, int major, int minor, int nfft=512, int nOverlap = 450, int padTo=800, int window=2, int log=0, int norm=0);

  TCanvas * avgRealNull(int ch, int major, int minor, int nfft=512, int nOverlap = 450, int padTo=800, int window=2, int log=0, int norm=0, int cursorType=0);

  TCanvas* avgRealFull(int ch, int major, int minor, int nfft=512, int nOverlap = 450, int padTo=800, int window=2, int log=0, int cursorType=0);
  
  int drawAvgRealNullWithGeom(int ch, int major, int minor, int nfft=512, int nOverlap = 450, int padTo=800, int window=2, int log=0, int norm=0);
  int drawAvgRealNullFull(int ch, int major, int minor, int nfft=512, int nOverlap=450, int padTo=800, int window=2, int log=0, int norm=0);

  int drawAvgRealFull(int ch, int major, int minor, int nfft=512, int nOverlap = 450, int padTo=800, int window=2, int log=0);

  int drawAvgRealFullWithGeom(int ch, int major, int minor, int nfft=512, int nOverlap = 450, int padTo=800, int window=2, int log=0);

  int drawAvg(int ch, int major, int minor);
  int drawAvgHilbert(int ch, int major, int minor, int datset=0);
  //utility functions

  TNtuple * integrateAllWithSideband(int major, int minor, int channel, int dataset, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, double sbxmin, double sbxmax, double sbymin, double sbymax, int norm);
  vector<vector<double>> sidebandSubtractAll(int major, int minor, int channel, int dataset, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, int norm);
  vector<vector<double>> sidebandSubtractAllXAxis(int major, int minor, int channel, int dataset, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, int norm);
  vector<vector<double>> sidebandSubtractAllYAxis(int major, int minor, int channel, int dataset, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, int norm);
  TNtuple * sidebandSubtractAllTuple(int major, int minor, int channel, int dataset, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, int norm);

  TH2D * nullSubtractedSpectrogram(int ch, int major, int minor, int event, int nfft, int nOverlap, int padTo, int window, int log);
  TH2D * avgNullSubtractedSpectrogram(int ch, int major, int minor, int nfft, int nOverlap, int padTo, int window, int log);

  double getSignalTOA(int ch, int major, int minor, int dataset=1);
  double getSignalDTOA(int ch, int major, int minor, int dataset);

  TH1F* bootstrapSidebandSubtract(int major, int minor, int channel, int dataset, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, int n, int nbins, int binlow, int binhigh, double constant=1.);

  TH1F* bootstrapSidebandSubtractXAxis(int major, int minor, int channel, int dataset, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, int n, int nbins, int binlow, int binhigh, double constant=1.);
  TH1F* bootstrapSidebandSubtractYAxis(int major, int minor, int channel, int dataset, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, int n, int nbins, int binlow, int binhigh, double constant=1.);
  int bootstrapSidebandSubtractXAxisVersusPower(TGraph * gr, int major, int minor, int channel, int dataset, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, int n);
  vector<double> getCharge(int major, int minor);
  int bootstrap2D(TGraph *gr, vector<double> vec0, vector<double>vec1);
  int bootstrap2D(TGraph *gr, vector<double> vec0, vector<vector<double>> vec1, int n);



  
}
#endif
