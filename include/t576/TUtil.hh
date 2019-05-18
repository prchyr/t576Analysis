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
  static constexpr double kB=8.617343e-11;//MeV/kelvin
  static constexpr double rho=1.168e-3;//sea level density
  static constexpr double x_0=36.7;//radiation length in air
  static constexpr double e_0=.078;//ionization energy 
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
  //approximated (fast) sinc interpolation.(will not be very good for low N.) sacrifice accuracy for speed.
  TGraph * sincInterpolateGraphFast(TGraph *inGr, double interpGSs, int N=10);
  /*interpolate a tgraph.
    
    type 0 is a ROOT akima interpolation.(default)
    type 1 is a sinc interpolation
    type 2 is a fast sinc interpolation. not as accurate as the sinc, but faster. the parameter N is only used for this type.

    the sinc interpolation is best/most accurate for truly band-limited signals sampled near the nyquist frequency. the fast method works very well, and is usually recommended over the full sinc for most applications. the akima is better for oversampled signals where the frequency of interest is far below the nyquist frequency.
  */
  int getInterpolatedGraph(TGraph * inGraph, TGraph *outGraph, double interpGSs, int type=0, int N=10);
  //return the interpolated graph (memory use) using ROOT's akima method
  TGraph * interpolateGraph(TGraph * inGraph, double interpGSs);

  //normalize a graph
  TGraph * normalize(TGraph *inGr);
  //normalize to peak
  TGraph * normToPeak(TGraph *inGr);
  //normalize a 2d graph
  TGraph2D * normalize(TGraph2D * inGr);
  //return a chunk of a graph, specified by x-axis values.
  //if shift_to_zero==1, the time axis is shifted such that it starts at 0.
  //  return the cumulative distribution function of a graph. if normed is 1, the graph is scaled such that x and y axes spread from 0 to 1.
  TGraph *CDF(TGraph *inGr, int normed=0);
  //get a chunk of a graph from start to end. delay to zero shifts the time so that the returned graph starts at t=0
  //uses TGraph->Eval() to get a very specific point in the graph. slow.
  TGraph * getChunkOfGraph(TGraph *ingr, double start, double end, int delay_to_zero=0);
//get a chunk of a graph from start to end. delay to zero shifts the time so that the returned graph starts at t=0
  TGraph * getChunkOfGraphFast(TGraph *ingr, double start, double end, int delay_to_zero=0);
  TGraph * getNSamplesFrom(TGraph *ingr, double start, int nSamples, int delay_to_zero);
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

//return the average of a bunch of aligned graphs
  TGraph* alignMultipleAndAverage(vector<TGraph*> inGr, double max_delay=999999., double t_low=0., double t_high=999999.);
  //align a large number of graphs to the first, then truncate them all to t_min and t_max
  vector<TGraph*> alignMultipleAndTruncate(vector<TGraph*> inGr, double max_delay, double t_min, double t_max, double t_low=0., double t_high=999999.);
  //align a large set to a reference set.
  vector<TGraph*> alignMultipleToOther(vector<TGraph*> inGr, vector<TGraph*> othGr, double max_delay=999999., double t_low=0., double t_high=999999.);
  //delay a graph
  TGraph *delayGraph(TGraph *ingr, double delay);
  //same but with no mem usage.
  int delayGraph(TGraph * ingr, TGraph *outgr, double delay);
  //plot \Delta(gr1[i], gr2[i]) for each graph point i
  TH1F * plotResiduals(TGraph *gr1, TGraph *gr2, int nbins=40, double min=1, double max=-1);
  //average some graphs
  TGraph * avgGraph(vector<TGraph*> inGr);
  //get the absolute value of a graph
  TGraph * absGraph(TGraph *inGr);
  //add 2 TGraphs. if constant is -1, they are subtracted.
  TGraph * add(TGraph * g1, TGraph * g2, double constant=1.);
  TGraph2D * add(TGraph2D * g1, TGraph2D * g2, double constant=1.);
  //dot product of 2 graphs
  double dot(TGraph *g1, TGraph *g2);
  //multiply two graphs: out(t)=g1(t)*consant*g2(t)
  TGraph * mult(TGraph *g1, TGraph *g2, double constant=1.);
  //divide 2 graphs: out(t)=g1(t)/constant*g2(t). if there is a divide
  //by zero, that entry will just be 0.
  TGraph * divide(TGraph *g1, TGraph *g2, double constant);
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
  //get the max value (wrapper of TMath::max)
  double max(TGraph *gr);
  double maxInRange(TGraph *gr, double t_low=0., double t_high=999999.);
  //get the x-axis location of t max value (wraper of TMath::LocMax)
  double locMax(TGraph *gr);
  double locMaxInRange(TGraph *gr, double t_low=0., double t_high=999999.);
  //get the index of a certain value. will always under-estimate
  int getIndex(TGraph *gr, double t);
  //get the power (square all values of a graph)
  TGraph * power(TGraph *gr);
  //remove the mean of a TGraph. range to compute the mean over is optional.
  //the mean computed within a sub-range will be removed from the full graph.
  TGraph * removeMean(TGraph *gr, double t_low=0., double t_high=999999.);
  //flip an array ([0]=[n-1]) (NOT WORKING)
  double * flip(int n, double * in);
  //flip a graph
  TGraph * flip(TGraph *inGr);
  //swap the x and y arrays of a TGraph
  TGraph * swap(TGraph * inGr);
  //same as removeMean but the original graph is changed
  int removeMeanInPlace(TGraph *gr, double t_low=0., double t_high=999999.);
  //make CW with given parameters.
  TGraph * makeCW(double freq,  double amp, double t_min=0., double t_max=1000., double GSs=20., double phase=0.);
  //integrate a TGraph. lower and upper bounds are optional.
  double integrate(TGraph * gr, double t_low=0, double t_high=999999.);  
  //get the RMS
  double rms(TGraph * gr, double t_low, double t_high);
//take the derivative. if direction=-1, takes derivative along other direction of axis.
  TGraph * derivative(TGraph *gr, int direction=1);
  //get the observer graph from a retarded graph. tUnits is the multiplier for ns (eg for ms, tUinits=1e6);
  TGraph * gObs(TGraph *inGr, double thetaDeg, double tUnits=1.);
  //integrate a TGraph but square it first to get units of power.
  double integratePower(TGraph * gr, double t_low=0, double t_high=999999.);
  //integrate a histogram
  double integrate(TH2D *h, double xmin, double xmax, double ymin, double ymax, double & err);
  double integrate2D(TH2D *h, double xmin, double xmax, double ymin, double ymax, double & err);
  //integrate a graph bin-by-bin, putting the result in a tgraph
  //binNS is the desired bin width in nanoseconds
  TGraph * integrateByBin(TGraph *gr, double binNS);
  //simple 2 pole lowpass filter
  TGraph * lowpassFilter(TGraph *ingr, double cutoff, int order=2);
  //highpass filter
  TGraph * highpassFilter(TGraph *ingr, double cutoff, int order=1);
  //a bandpass filter made of one of each of the above filters.
  TGraph * bandpassFilter(TGraph *ingr, double low, double high);
  //a brick wall frequency domain filter
  TGraph * brickWallFilter(TGraph *ingr, double low, double high);
    //will take the first chunk of the signal graph (equal to to t_high-t_low)

  // return the value of a window over sample numbers n. types are:
  /*
    0=bartlett (triangle) (default)
    1=welch (parabolic)
    2=hann (gaussian ish)
    3=blackman-nuttall (gaussian ish);
   */
  double window(int i, int n, int type=0);
  //return the value of a bartlett window over sample numbers n
  double bartlettWindow(int i, int n);
  //return the value of a welch window over sample numbers n
  double welchWindow(int i, int n);
  //return the value of a hann window over sample numbers n
  double hannWindow(int i, int n);
  //return the value of a blackman-nuttall window
  double blackmanNuttallWindow(int i, int n);
  //apply a window of the selected type to the graph inGr in time window.
  TGraph * applyWindow(TGraph* inGr, double startt, double endt, int type=0);
    //and add it to the indicated region of the background graph.
  TGraph * makeNullData(TGraph *sig, TGraph * back, double t_min, double t_max, double scale=1.);
  double sidebandSubtraction2DWithErrors(TH2D *h, double sband_x1, double sband_x2, double sband_y1, double sband_y2, double & err, int draw=0);
  double sidebandSubtraction2D(TH2D *h, double sband_x1, double sband_x2, double sband_y1, double sband_y2, int draw=0);
  //degrees to radians
  double deg2Rad(double deg);
  //radians to degrees
  double rad2Deg(double rad);
  //add some noise to a graph
  TGraph * addNoise(TGraph * inGr, double level);
  //find the x values of zero crossings. can also plot the relative time between subsequent zero crossings if relative = 1
TGraph * getZeroCrossGraph(TGraph * inGr, int relative=0);
  //get the time of the first threshold crossing after after
  double getFirstThresholdCrossing(TGraph *inGr, double thresh, double after=0.);


  //drawing things
 
  //draw a bunch of graphs
  void draw(vector<TGraph*> inGr, TString option="");
    //draw a bunch of graphs
  void draw(int nGraphs, TGraph ** inGr, TString option="");

  //a pretty warm palette
  void setWarmPalette();
  //a pretty cool palette
  void setCoolPalette();
  //cold palette
  void setColdPalette();
  void setHotPalette();
  void set2DPalette();
  //miscellaneous TGraph things


  //return a TGraph, evenly sampled along dt
  TGraph * evenSample(TGraph *inGr, double dt);
  //zero pad a tgraph. requires an evenly sampled tgraph.
  TGraph * zeroPad(TGraph *inGr, int num, int whichEnd=1);
  //get the distance between tvectors
  double distance3(TVector3 one, TVector3 two);
  //get the distance from one vector to several others
  double * distance3(int N, TVector3 one, TVector3 * two);
  //time of flight between two tvectors
  double timeOfFlight(TVector3 one, TVector3 two, double n=1.);
  //time of flight from one vector to several others
  double * timeOfFlight(int N, TVector3 one, TVector3 * two, double n=1.);
  //delta t between 2 antennas for a single source. it is the time a signal hits number 1 minus the time it hits number 2;
  double  dTimeOfFlight(TVector3 source,TVector3 one, TVector3 two, double n=1.);
  //delta t's between different antennas for a single source
  double ** dTimeOfFlight(int N, TVector3 one, TVector3 * two, double n=1.);
  //offset a bunch of antennas for the correct delta t's for single source
  //  TGraph ** delayGraphs
  namespace FFT{

    //returns a tgraph2d, x axis is the freqs, y axis is the real part of the complex fft, z axis is the imaginary part of the fft.
    TGraph2D * fft(TGraph *inGr);
    //returns the inverse fft of the data. must be structured as above (x axis is the freqs, y axis is the real part of the complex fft, z axis is the imaginary part of the fft.)
    TGraph * ifft(TGraph2D *inGr);
    //return the Hilbert transform
    TGraph * hilbertTransform(TGraph *inGr);
    //return the phase at a single point
    double phase(TGraph *inGr, double t);
    //plot the phase of the full graph (NOT WORKING)
    TGraph * plotPhase(TGraph *inGr);
    //zero the phase at the given frequency
    TGraph * zeroPhaseAt(TGraph *inGr, double freq, int debug=0);
    //set the phase angle at the given frequency
    TGraph * setPhaseAt(TGraph *inGr, double freq, double phaseAng, int debug=0);
    //return the Hilbert envelope
    TGraph * hilbertEnvelope(TGraph *inGr);
    //return the power spectral density in dBm/Hz, rBW is the resolution bandwith of the system, used to calculate the density. defaults to Nyquist.

    TGraph * psd(TGraph *inGr, double rBW=0., int dbFlag=1);


/*
      return the spectrogram, with various options:

binsize is what python calls nfft. it's the number of samples in the chunk over which the FFT is calculated. this is one spectrogram 'bin'

overlap is how many samples each bin overlaps with the next. helps somewhat with smoothing.

zero pad length is the length to which the chunk is symmetrically zero-padded.

win_type is an enumeration of window types to be applied to each bin. this helps avoid discontinuities and noise in the spectrogram. see the window function for the window types.
     */
    TH2D* spectrogram(TGraph *gr, Int_t binsize = 128, Int_t overlap=32, Int_t zero_pad_length=128, int win_type=0, int dbFlag=1);
    //averages a vector of spectrograms. must be the same size.
    TH2D* avgSpectrograms(vector<TH2D*>  inh);

//plot the peak frequency for each bin, with the options as above

TGraph* peakFreqGraph(TGraph *gr, Int_t binsize , Int_t overlap, Int_t zero_pad_length, int win_type, double thresh=0.);
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

    //must be for square matrix, a Ralston-style filter matrix
    TMatrixD makeFilter(TDecompSVD svd, int below, int above=0);
    //flatten the indices of a matrix in row major.
    TVectorD flatten(TMatrixD m);
    //cast a vector to graph for plotting.
    TGraph * toGraph(TVectorD v, double samplerate=1., double delay=0., TString name="");
    //makes a 2d histo from a matrix, for visualization.
    TH2F * matrixMap(TMatrixD M, TString name="");


  }

  namespace SIM{
    //the shower age expression for the NKG approximation  from arXiv:1503.02808
    double ss(double x, double E, double x_0, double e_0);
    //number of leptons as a function of depth in radiation lengths in the NKG approximation from  arXiv:1503.02808
    double n(double x, double E, double x_0, double e_0);
  }

}



#endif
