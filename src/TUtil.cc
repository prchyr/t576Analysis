/*
copyright s. prohira and the T576 collaboration 2018
released under the GNU General Public License version 3
*/

#include "TUtil.hh"

int TUtil::fN=1;
int TUtil::fNi=1;
TVirtualFFT * TUtil::fftr2c=TVirtualFFT::FFT(1, &fN, "R2C ES");
TVirtualFFT * TUtil::fftc2r=TVirtualFFT::FFT(1, &fN, "C2R ES");
TGraph2D *TUtil::fXfrmGr2D = new TGraph2D();
TGraph *TUtil::fXfrmGr = new TGraph();


/***************FFT things************************

the functions here utilise the above static members, meaning that
this class is very lightweight, and takes care of its own memory 
management. FFTW calculates the fastest and most optimal way to do
an FFT the first time it is initialized with new parameters. after that,
each subsequent fft of the same size is very fast. the outputs are stored
in the tgraphs above, so DO NOT delete the returned graph. 




**************************************************/

/*

the main FFT. returns a TGraph2D with x-axis=frequency, y-axis=real, 
and z-axis=imaginary. 

*/
TGraph2D * TUtil::FFT::fft(TGraph * inGr){
  int n = inGr->GetN();
  if(n!=fN){
    fN=n;
    fftr2c=TVirtualFFT::FFT(1, &fN, "R2C P K");
    //    fftc2r=TVirtualFFT::FFT(1, &fN, "C2R P");
  }
  
  double dt = inGr->GetX()[50]-inGr->GetX()[49];
  double fs = 1./dt;
  
  fftr2c->SetPoints(inGr->GetY());
  fftr2c->Transform();
  
  double re[n], im[n];
  fftr2c->GetPointsComplex(re, im);

  double df=fs/(n);

  double norm=1./sqrt((double)n);
  TGraph2D *outGr=new TGraph2D(n, makeIndices(n, df), &re[0], &im[0]);
  for (int i=0;i<outGr->GetN();i++){
    outGr->GetY()[i] *= norm;
    outGr->GetZ()[i] *= norm;
  }
  *fXfrmGr2D=*outGr;
  delete outGr;
  return fXfrmGr2D;
}



/*

Provide a TGraph 2D to this function and it will return the inverse fft
in the form of a tgraph.

 */


TGraph * TUtil::FFT::ifft(TGraph2D * inGr){
  int n = inGr->GetN();
  if(n!=fNi){
    fNi=n;
    //    fftr2c=TVirtualFFT::FFT(1, &fN, "R2C P");
    fftc2r=TVirtualFFT::FFT(1, &fNi, "C2R P K");
  }
  
  double df = inGr->GetX()[50]-inGr->GetX()[49];
  double fs=(double)n*df;
  double dt = 1./fs;
  
  fftc2r->SetPointsComplex(inGr->GetY(), inGr->GetZ());
  fftc2r->Transform();
  
  double re[n];
  fftc2r->GetPoints(re);
  double norm=1./sqrt((double)n);


  TGraph *outGr=new TGraph(n, makeIndices(n, dt), re);

  for (int i=0;i<outGr->GetN();i++) outGr->GetY()[i] *= norm;
  *fXfrmGr=*outGr;
  delete outGr;
  return fXfrmGr;
}

/*

psd, returns a tgraph in dbm/hz from an input tgraph

*/

TGraph * TUtil::FFT::psd(TGraph * inGr){
  auto xfrm=fft(inGr);
  int n=xfrm->GetN();
  //  double norm=1./(double)sqrt(n);
  auto xx=xfrm->GetX();
  auto re=xfrm->GetY();
  auto im=xfrm->GetZ();
  double yy[n];

  //resolution bandwidth, set to nyquist. should make this changeable. 
  double rbw=xx[n-1]/2;

  yy[0]=vToDbmHz(rbw,re[0]);
  for(int i=1;i<(n+1)/2;i++){
    yy[i]=vToDbmHz(rbw, re[i], im[i]);
  }
  yy[n/2]=vToDbmHz(rbw,re[n/2], im[n/2]);

  TGraph * outGr=new TGraph((n/2), xx, yy);

  *fXfrmGr=*outGr;
  delete outGr;
  return fXfrmGr;
}



TGraph * TUtil::FFT::hilbertTransform(TGraph *inGr){
  auto infft=fft(inGr);
  int n=infft->GetN();

  for(int i=0;i<n;i++){
    double im=infft->GetZ()[i];
    infft->GetZ()[i]=infft->GetY()[i];
    infft->GetY()[i]=-1.*im;
  }
  auto outGr=ifft(infft);
  return outGr;
}

TGraph * TUtil::FFT::hilbertEnvelope(TGraph * inGr){
  auto hilb=hilbertTransform(inGr);
  TGraph * out=new TGraph();
  for(int i=0;i<hilb->GetN();i++){
    out->SetPoint(i, inGr->GetX()[i], sqrt((inGr->GetY()[i]*inGr->GetY()[i])+(hilb->GetY()[i]*hilb->GetY()[i])));
  }
  *fXfrmGr=*out;
  delete out;
  return fXfrmGr;
}




/******************utility things*********************

none of these manage their own memory. normaliz() for example, will
make a new tgraph that you'll need to delete on your own.

*/

// TUtilGraph * TUtilGraph::operator*(const double a){
//   for(int i=0;i<this->GetN();i++)this->GetY()[i]*=a;
// }

double * TUtil::makeIndices(int n, double step, double offset){
  double *out=new double[n];
  for(int i=0;i<n;i++){
    out[i]=(i*step+offset);
  }
  return out;
}


double TUtil::vToDbmHz(double bandwidthGsNs, double re, double im){
  double val=re*re+im*im;
  return (10.*log10(val/50.))+30-(10.*log10(bandwidthGsNs*1.e9));
}

TGraph * TUtil::normalize(TGraph * inGr){
    double length=inGr->GetN();
    double norm=0;
    for(int i=0;i<length;i++){
      norm+=inGr->GetY()[i]*inGr->GetY()[i];
    }
    TGraph *og = (TGraph*)inGr->Clone();
    for(int i=0;i<length;i++)og->GetY()[i]*=norm;
    return og;
  }
  

int TUtil::getInterpolatedGraph(TGraph * inGraph, TGraph *outGraph, double interpGsNs){
  ROOT::Math::Interpolator interp(inGraph->GetN(), ROOT::Math::Interpolation::kAKIMA);
  interp.SetData(inGraph->GetN(), inGraph->GetX(), inGraph->GetY());

  //get dt, assuming even sampling.
  double inDt=inGraph->GetX()[50]-inGraph->GetX()[49];
  double inGsPerNs=1./inDt;

  double outDt=1./interpGsNs;

  int samps=(int) (inGraph->GetN()*(interpGsNs/inGsPerNs));

  vector<double> xx, yy;
  for(int i=0;i<samps;i++){
    double time = i*outDt;
    if(time>inGraph->GetX()[inGraph->GetN()-1])continue;
    xx.push_back(time);
    yy.push_back(interp.Eval(time));
  }

  auto tempGr=new TGraph(xx.size(), &xx[0], &yy[0]);
  *outGraph=*tempGr;
  delete(tempGr);
			 

  return 1;
}

TGraph * TUtil::getChunkOfGraphFine(TGraph *ingr, double start, double end){
  double *xx=ingr->GetX();
  double *yy=ingr->GetY();
  vector<double> outx, outy;
  double xincr=xx[10]-xx[9];
  for(int i=0;i<ingr->GetN();i++){
    //if(xx[i]>=start&&xx[i]<=end){
      //    }
      //else{
      double time=start+((double)i*xincr);
      if(time<end){
      outx.push_back(time);
      //      outx.push_back(xx[i]);
      outy.push_back(ingr->Eval(time));
      }
  }
  TGraph * outg=new TGraph(outx.size(), &outx[0], &outy[0]);
  return outg;
}
