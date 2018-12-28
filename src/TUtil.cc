#include "TUtil.hh"

int TUtil::fN=1;
TVirtualFFT * TUtil::fftr2c=TVirtualFFT::FFT(1, &fN, "R2C ES K");

TGraph * TUtil::FFT::psd(TGraph * inGr){
  int n = inGr->GetN();
  cout<<n<<endl;
  if(n!=fN){
    fN=n;
    fftr2c=TVirtualFFT::FFT(1, &fN, "R2C P K");
  }
  cout<<fN<<endl;
  double dt = inGr->GetX()[50]-inGr->GetX()[49];
  double fs = 1./dt;
  
  fftr2c->SetPoints(inGr->GetY());
  fftr2c->Transform();
  double re[n], im[n];
  fftr2c->GetPoints(re, im);

  double yy[n/2], xx[n/2];
  double df=fs/(2*n);
  for(int i=0;i<(n/2)-1;i++){
    xx[i]=i*df;
    yy[i]=sqrt(pow(re[i], 2)+pow(im[i], 2));
  }
  TGraph * outGr=new TGraph((n/2)-1, xx, yy);
  
  return outGr;
}
