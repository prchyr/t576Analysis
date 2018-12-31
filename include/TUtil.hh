/*
copyright s. prohira and the T576 collaboration 2018
released under the GNU General Public License version 3
*/

#ifndef TUTIL_BASE_H
#define TUTIL_BASE_H

#include "T576Event.hh"
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
#include "TDecompSVD.h"
#include <iostream>
#include <fstream>
#include <vector>

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "TUtilGraph.hh"


using namespace CLHEP;
using namespace std;

class TUtilGraph;
class TUtil : public TObject{
public:

  static double setN(double val){fN=val;return fN;};
  static double getN(){return fN;}
  static double vToDbmHz(double bandwidthGsNs, double re, double im=0);
  static double * makeIndices(int n, double step, double offset=0);
  static int getInterpolatedGraph(TGraph * inGraph, TGraph *outGraph, double interpGsNs);
  static TGraph * normalize(TGraph *inGr);
  static TGraph * getChunkOfGraphFine(TGraph *ingr, double start, double end);
  static TGraph *delayGraph(TGraph *ingr, double delay);
  static int delayGraph(TGraph * ingr, TGraph *outgr, double delay);
  static TH1F * plotResiduals(TGraph *gr1, TGraph *gr2, int nbins=40, double min=1, double max=-1);

  static void setWarmPalette();
  static void setCoolPalette();
  
  class FFT{
  public:
    //returns a tgraph2d, x axis is the freqs, y axis is the real part of the complex fft, z axis is the imaginary part of the fft.
    static TGraph2D * fft(TGraph *inGr);
    //returns the inverse fft of the data. must be structured as above (x axis is the freqs, y axis is the real part of the complex fft, z axis is the imaginary part of the fft.)
    static TGraph * ifft(TGraph2D *inGr);
    //return the Hilbert transform
    static TGraph * hilbertTransform(TGraph *inGr);
    //return the Hilbert envelope
    static TGraph * hilbertEnvelope(TGraph *inGr);
    //return the power spectral density in dBm/Hz
    static TGraph * psd(TGraph *inGr);
    //return the spectrogram, with various options. 
    static TGraph2D * spectrogram(TGraph *inGr, int nfft=256, int noverlap=0);

  private:

    ClassDefNV(FFT, 1);
  };

  class SVD{
  public:
    static TVectorD normalize(TVectorD vec);
    static double norm(TVectorD vec);
    //build a matrix out of events, with each event a row in the matrix.
    //must all be the same length.
    static TMatrixD eventMatrix(vector<TGraph*> vecs);
    //make a density matrix partitioned along D
    static TMatrixD densityMatrix(TGraph *vec, int D=1, int xlim=0);
    //returns a matrix after removing all coffs higher than below.
    //if above is provided, will return the matrix as constructed from
    //the modes between above and below.
    static TMatrixD truncateSVD(TDecompSVD svd, int below, int above=0);
    //returns a matrix constructed from a single singular value
    static TMatrixD reconstructSingle(TDecompSVD svd, int val);
    //cast a tgraph to a tvectord
    static TVectorD toVector(TGraph *vec);
    //makes a row average of matrix m
    static TVectorD avgVector(TMatrixD m);
    //builds a basis of patterns from the matrix m, up to the number of
    //patterns num
    static TMatrixD buildBasis(TMatrixD m, int num=10);
    //returns the expansion coeffs for a vector in a basis B
    static TVectorD getCoefficients(TVectorD V, TMatrixD B);
    //expands a vector in a basis B
    static TVectorD expandInBasis(TVectorD V, TMatrixD B, int num=10);
    static TVectorD expandInBasis(TGraph * G, TMatrixD B, int num=10);
    //must be for square matrix, a Ralston-style filter matrix
    static TMatrixD makeFilter(TDecompSVD svd, int below, int above=0);
    //flatten the indices of a matrix in row major.
    static TVectorD flatten(TMatrixD m);
    //cast a vector to graph for plotting.
    static TGraph * toGraph(TVectorD v, double samplerate=1., double delay=0., TString name="");
    //makes a 2d histo from a matrix, for visualization.
    static TH2F * matrixMap(TMatrixD M, TString name="");

  private:
    ClassDefNV(SVD, 1);
  };
  
private:

  static int fN;
  static int fNi;
  static TVirtualFFT *fftr2c;
  static TVirtualFFT *fftc2r;
  static TGraph2D *fXfrmGr2D;
  static TGraph * fXfrmGr;
  ClassDefNV(TUtil, 1);
};



#endif
