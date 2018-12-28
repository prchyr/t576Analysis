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


using namespace CLHEP;
using namespace std;


class TUtil : public TObject{
public:

  static double setN(double val){fN=val;return fN;};
  static double getN(){return fN;}
  static double vToDbmHz(double sampRateGsNs, double re, double im=0);
  
  class FFT{
  public:
    //returns a tgraph2d, x axis is the freqs, y axis is the real part of the complex fft, z axis is the imaginary part of the fft.
    static TGraph2D * fft(TGraph *inGr);
    //returns the inverse fft of the data. must be structured as above (x axis is the freqs, y axis is the real part of the complex fft, z axis is the imaginary part of the fft.)
    static TGraph * ifft(TGraph2D *inGr);
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
  static TVirtualFFT *fftr2c;
  static TVirtualFFT *fftc2r;
  ClassDefNV(TUtil, 1);
};



#endif
