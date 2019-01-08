/*
  copyright s. prohira and the T576 collaboration 2018
  released under the GNU General Public License version 3
*/

#ifndef TUTIL_BASE_H
#define TUTIL_BASE_H


#include "TROOT.h"
#include "TRint.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TObject.h"
#include "TFile.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "cnpy.h"
#include "TSystemFile.h"
#include "TH1F.h"
#include "TH2F.h"
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
#include "TMatrixD.h"
#include "TVector3.h"
#include "TDecompSVD.h"
#include <iostream>
#include <fstream>
#include <vector>

// #include "CLHEP/Units/PhysicalConstants.h"
// #include "CLHEP/Vector/LorentzVector.h"
// #include "CLHEP/Vector/ThreeVector.h"

#include "TUtilGraph.hh"

#define c_light .29979246 //  m/ns
#define pi 3.1415927

//using namespace CLHEP;
using namespace std;

class TUtilGraph;
namespace TUtil{

  //volts to dbm/hz
  double vToDbmHz(double bandwidthGSs, double re, double im=0);
  //make an axis with linearly increasing values.
  double * makeIndices(int n, double step, double offset=0);
  //the normalized sinc function: sin(pi x)/(pi x)
  double sinc(double x);
  //return a TGraph interpoltaed using simple sinc interpolation.
  TGraph * sincInterpolateGraph(TGraph *inGr, double interpGSs);
  TGraph * sincInterpolateGraphSimple(TGraph *inGr, double interpGSs);
  //interpolate a tgraph using the ROOT interpolation functions
  int getInterpolatedGraph(TGraph * inGraph, TGraph *outGraph, double interpGSs);

  //normalize a graph
  TGraph * normalize(TGraph *inGr);
  //return a chunk of a graph, specified by x-axis values. 
  TGraph * getChunkOfGraph(TGraph *ingr, double start, double end);
  //cross correlation of two graphs. returns the cross-correlation graph
  //maxDelay is the maximum starting offset between gr1 and gr2. defaults
  //to the full length of the graphs.
  //t_low(high) are the lowest(highest) times over which to
  //calculate the correlation. defaults to the full length of the graphs.
  TGraph * crossCorrelate(TGraph * gr1, TGraph * gr2, double max_delay=999999., double t_low=0., double t_high=999999.);
  //same as above but allows you to supply a window function, to only
  //correlate parts of the graph that you want.
  TGraph * crossCorrelateWindowed(TGraph * gr1, TGraph * gr2, TGraph *wingraph, double max_delay=999999., double t_low=0., double t_high=999999.);
  //the same as the crossCorrelate() function, but returns gr2 shifted in time
  //to the point of peak cross correlation with gr1.
  TGraph * align(TGraph * gr1, TGraph * gr2, double max_delay=999999., double t_low=0., double t_high=999999.);
  //align gr2 to gr1, but returning othGr, delayed correctly. 
  TGraph * alignToOther(TGraph * gr1, TGraph * gr2, TGraph* othGr, double max_delay=999999., double t_low=0., double t_high=999999.);
  //align a large number of graphs to the first graph in the set.
  vector<TGraph*> alignMultiple(vector<TGraph*> inGr, double max_delay=999999., double t_low=0., double t_high=999999.);
  //align a large set to a reference set.
  vector<TGraph*> alignMultipleToOther(vector<TGraph*> inGr, vector<TGraph*> othGr, double max_delay=999999., double t_low=0., double t_high=999999.);
  //delay a graph
  TGraph *delayGraph(TGraph *ingr, double delay);
  //same but with no mem usage.
  int delayGraph(TGraph * ingr, TGraph *outgr, double delay);
  //plot \Delta(gr1[i], gr2[i]) for each graph point i
  TH1F * plotResiduals(TGraph *gr1, TGraph *gr2, int nbins=40, double min=1, double max=-1);
  //add 2 TGraphs. if constant is -1, they are subtracted.
  TGraph * add(TGraph * g1, TGraph * g2, double constant=1.);
  double integrate(TGraph * gr, double t_low=0, double t_high=999999.);
  double deg2Rad(double deg);
  double rad2Deg(double rad);
  void setWarmPalette();
  void setCoolPalette();
  
  namespace FFT{

    //returns a tgraph2d, x axis is the freqs, y axis is the real part of the complex fft, z axis is the imaginary part of the fft.
    TGraph2D * fft(TGraph *inGr);
    //returns the inverse fft of the data. must be structured as above (x axis is the freqs, y axis is the real part of the complex fft, z axis is the imaginary part of the fft.)
    TGraph * ifft(TGraph2D *inGr);
    //return the Hilbert transform
    TGraph * hilbertTransform(TGraph *inGr);
    //return the Hilbert envelope
    TGraph * hilbertEnvelope(TGraph *inGr);
    //return the power spectral density in dBm/Hz
    TGraph * psd(TGraph *inGr);
    //return the spectrogram, with various options. 
    TGraph2D * spectrogram(TGraph *inGr, int nfft=256, int noverlap=0);


  }

  namespace SVD{

    //normmalize a tvector
    TVectorD normalize(TVectorD vec);
    double norm(TVectorD vec);
    //build a matrix out of events, with each event a row in the matrix.
    //must all be the same length.
    TMatrixD eventMatrix(vector<TGraph*> vecs);
    //make a density matrix partitioned along D
    TMatrixD densityMatrix(TGraph *vec, int D=1, int xlim=0);
    //returns a matrix after removing all coffs higher than below.
    //if above is provided, will return the matrix as constructed from
    //the modes between above and below.
    TMatrixD truncateSVD(TDecompSVD svd, int below, int above=0);
    //returns a matrix constructed from a single singular value
    TMatrixD reconstructSingle(TDecompSVD svd, int val);
    //cast a tgraph to a tvectord
    TVectorD toVector(TGraph *vec);
    //makes a row average of matrix m
    TVectorD avgVector(TMatrixD m);
    //builds a basis of patterns from the matrix m, up to the number of
    //patterns num
    TMatrixD buildBasis(TMatrixD m, int num=10);
    //returns the expansion coeffs for a vector in a basis B
    TVectorD getCoefficients(TVectorD V, TMatrixD B);
    //expands a vector in a basis B
    TVectorD expandInBasis(TVectorD V, TMatrixD B, int num=10);
    TGraph * expandInBasis(TGraph * G, TMatrixD B, int num=10);
    TVectorD filter(TVectorD V, TMatrixD B, int num);
    TGraph * filter(TGraph *G, TMatrixD B, int num);

    //must be for square matrix, a Ralston-style filter matrix
    TMatrixD makeFilter(TDecompSVD svd, int below, int above=0);
    //flatten the indices of a matrix in row major.
    TVectorD flatten(TMatrixD m);
    //cast a vector to graph for plotting.
    TGraph * toGraph(TVectorD v, double samplerate=1., double delay=0., TString name="");
    //makes a 2d histo from a matrix, for visualization.
    TH2F * matrixMap(TMatrixD M, TString name="");


  }
  

}



#endif
