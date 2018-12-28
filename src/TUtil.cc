#include "TUtil.hh"

int TUtil::fN=1;
TVirtualFFT * TUtil::fftr2c=TVirtualFFT::FFT(1, &fN, "R2C ES");



TGraph * TUtil::FFT::psd(TGraph * inGr){
  int n = inGr->GetN();
  cout<<n<<" "<<fN<<" "<<*fftr2c->GetTransformFlag()<<endl;
  if(n!=fN){
    fN=n;
    fftr2c=TVirtualFFT::FFT(1, &fN, "R2C P");
  }
  //  fftr2c->Init("EX");

  cout<<n<<" "<<fN<<" "<<*fftr2c->GetTransformFlag()<<endl;
  double dt = inGr->GetX()[50]-inGr->GetX()[49];
  double fs = 1./dt;
  cout<<dt<<" "<<fs<<" "<<endl;
  fftr2c->SetPoints(inGr->GetY());
  fftr2c->Transform();
  //  cout<<*fftr2c->GetN()<<endl;
  double re[n], im[n];
  fftr2c->GetPointsComplex(re, im);

  double yy[n/2], xx[n/2];
  double df=fs/(n);

  //resolution bandwidth
  double rbw=fs/2.;
  double norm=1./(double)sqrt(n);
  cout<<df<<endl;
  xx[0]=0;
  yy[0]=vToDbmHz(rbw,norm*re[0]);
  for(int i=1;i<(n+1)/2;i++){
    xx[i]=(double)i*df;
    yy[i]=vToDbmHz(rbw, norm*re[i], norm*im[i]);
  }
  xx[n/2]=(double)n*df/2.;
  yy[n/2]=vToDbmHz(rbw,norm*re[n/2], norm*im[n/2]);
  
  TGraph * outGr=new TGraph((n/2), xx, yy);
  
  return outGr;
}

double TUtil::vToDbmHz(double sampRateGsNs, double re, double im){
  double val=re*re+im*im;
  return (10.*log10(val/50.))+30-(10.*log10(sampRateGsNs*1.e9));
}
