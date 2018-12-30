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
#include "TString.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "TVirtualFFT.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"
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
