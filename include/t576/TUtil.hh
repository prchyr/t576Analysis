/*
  copyright s. prohira and the T576 collaboration 2018
  released under the GNU General Public License version 3
*/


/*
The TUtil namespace provides numerous useful functions which act primarily on TGraphs. 

The FFT namespace has FFT functions built on ROOT's implementation of FFTW3. 

The SVD namespace has SVD functions built on the ROOT implementation of the GSL linear algebra tools.

The global system of units is: 

time: ns
length: m
frequency: GHz
power spectral density: dBm/Hz
amplitude: V
charge: nC

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
#include "TMatrixD.h"
#include "TVector3.h"
#include "TDecompSVD.h"
#include <iostream>
#include <fstream>
#include <vector>



#include "TUtilGraph.hh"



using namespace std;




class TUtilGraph;
namespace TUtil{

  /*some useful units. inspired by the CLHEP global system of units. 

    use these to keep your numbers in the global system of units defined above:
    ns, GHz, m, nC

    if you want to write something in terms of MHz, for example, just do 
    
    freq = 1200*MHz

    and then freq will have units of GHz, as it should.
  */
  //constants
  static constexpr double c_light = .29979246; //  m/ns
  static constexpr double pi = 3.1415927; //radians
  static constexpr double z_0=50; //ohms
  static constexpr double deg=pi/180.; //radians
  /*lengths
 
    for example, if you wanted to calculate the time it took for a signal 
    to propagate 75 feet, you'd do:
    
    75*ft/c_light
    
    and it would return the correct time in nanoseconds. 
  */
  static constexpr double m = 1.;
  static constexpr double ft = .3047*m;
  static constexpr double cm = .01*m;
  static constexpr double mm = .001*m;


  //time
  static constexpr double ns = 1.;
  static constexpr double us = ns*1e3;
  static constexpr double ms = ns*1e6;
  static constexpr double s = ns*1e9;

  //frequency
  static constexpr double GHz = 1.;
  static constexpr double MHz = .001*GHz;
  static constexpr double kHz = 1e-6*GHz;
  static constexpr double Hz = 1e-9*GHz;
  

  
  //volts to dbm/hz
  double vToDbmHz(double bandwidthGSs, double re, double im=0);
  //make an axis with linearly increasing values.
  double * makeIndices(int n, double step, double offset=0);
  //the normalized sinc function: sin(pi x)/(pi x)
  double sinc(double x);
  //return a TGraph interpoltaed using simple sinc interpolation.(slow)
  TGraph * sincInterpolateGraph(TGraph *inGr, double interpGSs);
  //approximated (fast) sinc interpolation.(broken)
  TGraph * sincInterpolateGraphSimple(TGraph *inGr, double interpGSs);
  //interpolate a tgraph using the ROOT interpolation functions(works well)
  int getInterpolatedGraph(TGraph * inGraph, TGraph *outGraph, double interpGSs);

  //normalize a graph
  TGraph * normalize(TGraph *inGr);
  //return a chunk of a graph, specified by x-axis values.
  //if shift_to_zero==1, the time axis is shifted such that it starts at 0.
  TGraph * getChunkOfGraph(TGraph *ingr, double start, double end, int delay_to_zero=0);
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
  /*align gr2 to gr1, but returning othGr, delayed by the amount needed to
  align gr2 to gr1 at the point of peak correlation. 
  for example, gr1 and gr2 are events which are causal, and othGr
  is some other graph which you'd like to align with these, but can't for
  whatever reason (contaminated with CW, etc.). 
  */
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
  //shift a graph along the y axis by the factor
  TGraph * shiftY(TGraph *g1, double factor);
  //shift a graph along the x axis by the factor
  TGraph * shiftX(TGraph *g1, double factor);
  //scale a TGraph by a constant factor
  TGraph * scale(TGraph *g1, double factor);
  //stretch a TGraph in time by a factor
  TGraph *stretch(TGraph *g1, double factor);
  //find the mean of a TGraph. range is optional
  double mean(TGraph *gr, double t_low=0., double t_high=999999.);
  //remove the mean of a TGraph. range to compute the mean over is optional.
  //the mean computed within a sub-range will be removed from the full graph.
  TGraph * removeMean(TGraph *gr, double t_low=0., double t_high=999999.);
  //make CW with given parameters.
  TGraph * makeCW(double freq,  double amp, double t_min=0., double t_max=1000., double GSs=20., double phase=0.);
  //integrate a TGraph. lower and upper bounds are optional.
  double integrate(TGraph * gr, double t_low=0, double t_high=999999.);
  //simple 2 pole lowpass filter
  TGraph * lowpassFilter(TGraph *ingr, double cutoff, int order=2);
  //degrees to radians
  double deg2Rad(double deg);
  //radians to degrees
  double rad2Deg(double rad);
  //a pretty warm palette
  void setWarmPalette();
  //a pretty cool palette
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
    //return the power spectral density in dBm/Hz, rBW is the resolution bandwith of the system, used to calculate the density. defaults to Nyquist.
    TGraph * psd(TGraph *inGr, double rBW=0.);
    //return the spectrogram, with various options. 
    TH2D* spectrogram(TGraph *gr, Int_t binsize = 128, Int_t overlap=32, Int_t zero_pad_length=128);
    //averages a vector of spectrograms. must be the same size.
    TH2D* avgSpectrograms(vector<TH2D*>  inh);

  }

  namespace SVD{

    //normmalize a tvector
    TVectorD normalize(TVectorD vec);
    //return the norm of a vector
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

    //filter a vector using a basis. the expansion of the vector in the basis B (to
    //order num) will be removed from the vector.
    TVectorD filter(TVectorD V, TMatrixD B, int num);
    TGraph * filter(TGraph *G, TMatrixD B, int num);
    //will align a signal matrix and a background matrix to their respective
    //reference matrices, build a basis out of backs ,
    //and filter sigs using this basis to order num
    TGraph * alignAndFilter(vector<TGraph*> sigs, vector<TGraph*> sigsref, vector<TGraph*> backs, vector<TGraph*> backsref, int num);
    //will take the first chunk of the signal graph (equal to to t_high-t_low)
    //and add it to the indicated region of the background graph.
    TGraph * makeNullData(TGraph *sig, TGraph * back, double t_min, double t_max);
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
