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

// TH2D* spectrogram(TGraph *gr, Int_t binsize = 128, Int_t overlap=32, Int_t zero_pad_length=128, int log=1, int draw=1){

//   Int_t size = gr->GetN();
//   double samprate=1./gr->GetX()[1]-gr->GetX()[0]);
//   double xmax = size/samprate;
//   double xmin = 0;

//   Int_t num_zeros=(zero_pad_length-binsize)/2.;
//   Int_t nbins = size/overlap;
//   char*timebuff;
//   double samplerate = size/(xmax-xmin);
//   double bandwidth = 1e9*samplerate;
//   TH1F *in = new TH1F("inhi", "inhi", zero_pad_length, 0, zero_pad_length);  
//   TH1*outt=0;
//   //TH2F *outhist=new TH2F("outhist", "spectrogram", nbins, xmin, xmax, (binsize), 0, samplerate);
//   //  cout<<binsize<<" "<<samplerate<<endl;
// TH2D *spectrogramHist=new TH2D("outhist", "spectrogram", nbins, xmin, xmax, (zero_pad_length), 0, samplerate);

//   Int_t start = 0;
//   //  Int_t j=0;
// TGraph *fGr=new TGraph(binsize+num_zeros+num_zeros);
//   for(int i=0;i<=nbins;i++){
//     for(int j=0;j<num_zeros;j++){
//     //    for(int j = 0;j<=zero_pad_length;j++){
//       in->SetPoint(j, 0);
//     }
//     for(int j=num_zeros;j<binsize+num_zeros;j++){
//       if((j+start)>=size)break;
//       in->SetBinContent(j, data[j+start]*FFTtools::bartlettWindow(j-num_zeros, binsize));
//     }
//     for(int j=binsize+num_zeros;j<zero_pad_length;j++){
//       in->SetBinContent(j, 0);
//     }

//     outt=in->FFT(outt, "mag");
//     outt->Scale(1./sqrt(zero_pad_length));
//     for(int j = 0;j<=(zero_pad_length);j++){
//       Double_t y = outt->GetBinContent(j);
//       if(log==1){
// 	y = (10.*log10(pow(y, 2.)/50.));//mv->dbm
// 	//y=10.*log10(y)+30;
// 	spectrogramHist->SetBinContent(i,j,(y-(10.*log10(bandwidth/2.))));//dmb/hz
//       }
//       else{
//       spectrogramHist->SetBinContent(i,j,y);//dmb
//       }
//     }
//     start+=overlap-1;
//   }
//   //cout<<"here"<<endl;
//   //  ->SetRightMargin(.15);
//   if(draw==1){
//   spectrogramHist->GetYaxis()->SetRangeUser(0, spectrogramHist->GetYaxis()->GetXmax()/2);
  
//   spectrogramHist->GetXaxis()->SetTitle("Time (ns)");

//   spectrogramHist->GetYaxis()->SetTitle("Frequency (GHz)");

//   spectrogramHist->GetZaxis()->SetTitle("dBm/Hz");
//   spectrogramHist->GetZaxis()->SetTitleOffset(1.5);
//   spectrogramHist->SetStats(0);
//   spectrogramHist->Draw("colz");
//   //maxbins->Draw("same");
//   //maxvals->SetLineColor(kRed);
//   //maxvals->Draw("same");
//   }
//   outt->Delete();
//   in->Delete();

//   // spectrogramHist->Delete();
//   return spectrogramHist;
//   //return outdat;
// }






/*******************SVD things**********************

these functions use the root linear algebra routines, which are 
quite extensive. they are built on the usual GSL routines. 

memory management is not done here.


 */

TVectorD TUtil::SVD::normalize(TVectorD vec){
  auto a=TVectorD(vec);
  auto b=TVectorD(vec);
  int len=vec.GetNrows();
  double norm=0;
  a.Sqr();
  norm=sqrt(a.Sum());
  b*=(1./norm);
  return b;
}


double TUtil::SVD::norm(TVectorD vec){
  auto a=TVectorD(vec);
  int len=vec.GetNrows();
  double norm=0;
  a.Sqr();
  norm=sqrt(a.Sum());
  return norm;
}

TMatrixD TUtil::SVD::eventMatrix(vector<TGraph*> vecs){
  int x=vecs.size();
  int y=vecs[0]->GetN();
  TMatrixD M(x,y);
  M.Zero();
  for(int i=0;i<x;i++){
    for(int j=0;j<y;j++){
    M[i][j]=vecs[i]->GetY()[j];
    }
  }
  return M;
}

//make a density matrix partitioned along D
TMatrixD TUtil::SVD::densityMatrix(TGraph *vec, int D, int xlim){
  int N=xlim==0?vec->GetN():xlim;
  int d=(int) (N/D);
  TVectorD V(vec->GetN(), vec->GetY());
  TMatrixD A(D, D);
  A.Zero();
  for(int i=0;i<d;i++){
    auto vv=V.GetSub(i*D, (i*D)+D, "");
    A.Rank1Update(vv, 1.);
  }
  return A;
}

TMatrixD TUtil::SVD::truncateSVD(TDecompSVD svd, int below, int above){
  auto V=svd.GetV();
  auto S=svd.GetSig();
  auto U=svd.GetU();
  auto temp=TMatrixD(U.GetNcols(), V.GetNcols());
  //cout<<U.GetNcols()<<" "<<V.GetNcols()<<endl;
  //cout<<temp.GetNrows()<<" "<<temp.GetNcols()<<endl;

  temp.Zero();
  for(int i=above;i<below;i++){
    temp[i][i]=S[i];
  }
  //cout<<"here"<<endl;
  auto SV=temp*V.T();
  //cout<<"no here"<<endl;
  //cout<<temp.GetNrows()<<" "<<temp.GetNcols()<<endl;
  auto USV=U*temp;
  //  cout<<"no no here"<<endl;
  return USV;
}


TMatrixD TUtil::SVD::reconstructSingle(TDecompSVD svd, int val){
  auto M = truncateSVD(svd, val+1, val);
  return M;
}

TVectorD TUtil::SVD::toVector(TGraph *vec){
  auto v=TVectorD(vec->GetN(), vec->GetY());
  return v;
}

//TVectorD subtract(

//TVectorD hilbertTransform(TVectorD vec){
  

TVectorD TUtil::SVD::avgVector(TMatrixD m){
  int sizex=m.GetNcols();
  int sizey=m.GetNrows();
  auto vec=TVectorD(sizex);
  //vector<double>vec;
  for(int i=0;i<sizey;i++){
    vec+=m[i];
  }
  return vec;
}  


TMatrixD TUtil::SVD::buildBasis(TMatrixD m, int num){
  auto rows=m.GetNrows();
  auto cols=m.GetNcols();
  //cout<<rows<<" "<<cols<<endl;
  num=num>cols?cols:num;
  auto M=TMatrixD(num, cols);

  auto mm=rows>=cols?m:m.T();
  //  cout<<mm.GetNrows()<<" "<<mm.GetNcols()<<endl;
  auto svd=TDecompSVD(mm);
  rows=M.GetNrows();
  cols=M.GetNcols();
  //cout<<rows<<" "<<cols<<endl;
  for(int i=0;i<rows;i++){
    auto recon=reconstructSingle(svd,i).T();
    //cout<<recon.GetNrows()<<" "<<recon.GetNcols()<<endl;
    auto avg=avgVector(recon);
    //cout<<avg.GetNrows()<<endl;//" "<<avg.GetNcols()<<endl;
    auto vec=normalize(avg);
    //cout<<vec.GetNrows()<<endl;
    for(int j=0;j<cols;j++){
      M[i][j]=vec[j];
    }
  }
  return rows>cols?M.T():M;
}

TVectorD TUtil::SVD::getCoefficients(TVectorD V, TMatrixD B){
  int y=B.GetNrows();
  int x=B.GetNcols();
  auto vec=TVectorD(y);
  vec.Zero();
  for (int i=0;i<y;i++){
    auto vnorm=normalize(V);
    for( int j=0;j<x;j++){
      vnorm[j]*=B[i][j];
    }
    vec[i]=vnorm.Sum();
  }
  return vec;
}

TVectorD TUtil::SVD::expandInBasis(TVectorD V, TMatrixD B, int num){

  int y=B.GetNrows();
  int x=B.GetNcols();
  auto vec=getCoefficients(V, B);
  auto outvec=TVectorD(V.GetNrows());
  outvec.Zero();
  num=num>=y?y:num;
  for (int i=0;i<num;i++){
    //    auto aligned=
    for (int j=0;j<x;j++){
      outvec[j]+=(vec[i]*B[i][j]);
    }
  }
  outvec*=norm(V);
  return outvec;
}

TVectorD TUtil::SVD::expandInBasis(TGraph * G, TMatrixD B, int num){
  TVectorD V = toVector(G);
  int y=B.GetNrows();
  int x=B.GetNcols();
  auto vec=getCoefficients(V, B);
  auto outvec=TVectorD(V.GetNrows());
  outvec.Zero();
  num=num>=y?y:num;
  for (int i=0;i<num;i++){
    //    auto aligned=
    for (int j=0;j<x;j++){
      outvec[j]+=(vec[i]*B[i][j]);
    }
  }
  outvec*=norm(V);
  return outvec;
}


//must be for square matrix
TMatrixD TUtil::SVD::makeFilter(TDecompSVD svd, int below, int above){
  auto M = truncateSVD(svd, below, above);
  auto D = M.GetNrows();
  auto I = TMatrixD(D, D);
  I.UnitMatrix();
  auto filter=I-(M*(1./M.NormInf()));
  //  auto filter=(M*(1./M.E2Norm()));
  return filter;
}


TGraph * TUtil::SVD::toGraph(TVectorD v, double samplerate, double delay,TString name){
  double tdiv=1./samplerate;
  auto  gg=new TGraph(v.GetNrows(), makeIndices(v.GetNrows(), tdiv, delay), &v[0]);
  gg->SetTitle("");
  gg->SetName(name);
  gg->GetXaxis()->SetRangeUser(0,gg->GetX()[gg->GetN()-1]);
  return gg;
}

TVectorD TUtil::SVD::flatten(TMatrixD m){
  int sizex=m.GetNcols();
  int sizey=m.GetNrows();
    auto vec=TVectorD(sizex*sizey);
  //vector<double>vec;
  for(int i=0;i<sizey;i++){
    vec.SetSub(i*sizex, m[i]);
  }
  return vec;
}  

TH2F * TUtil::SVD::matrixMap(TMatrixD M, TString name){
  int sizex=M.GetNcols();
  int sizey=M.GetNrows();
  //  cout<<sizey<<" "<<sizex<<endl;
  TH2F * map=new TH2F(name, name, sizex, 0.,sizex, sizey, 0,sizey);
  for(int i=0;i<sizex;i++){
    for(int j=0;j<sizey;j++){
      map->SetBinContent(i+1,j+1, M[i][sizey-j-1]);
    }
  }
  
  return map;
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



TGraph * TUtil::delayGraph(TGraph *ingr, double delay){
  double*xx=ingr->GetX();
  double xxout[ingr->GetN()];
  double*yy=ingr->GetY();
  for(int i=0;i<ingr->GetN();i++){
    xxout[i]=(double)xx[i]+delay;
  }
  TGraph *dg=new TGraph(ingr->GetN(), xxout, yy);
  return dg;
}

int TUtil::delayGraph(TGraph *ingr, TGraph *outgr, double delay){
  double*xx=ingr->GetX();
  double xxout[ingr->GetN()];
  double*yy=ingr->GetY();
  for(int i=0;i<ingr->GetN();i++){
    xxout[i]=(double)xx[i]+delay;
  }
  TGraph *dg=new TGraph(ingr->GetN(), xxout, yy);
  *outgr=*dg;
  delete dg;
  return 1;
}


TH1F * TUtil::plotResiduals(TGraph *gr1, TGraph *gr2, int nbins, double min, double max){
  TH1F *hist =new TH1F("", "", nbins,min, max);
  for(int i=0;i<gr1->GetN();i++)hist->Fill(gr1->GetY()[i]-gr2->GetY()[i]);
  return hist;
}


TGraph * crossCorrelate(TGraph * xx, TGraph * yy, int max_delay=100000, int xcorr_graph=0){
  double *x = xx->GetY();
  double *time=xx->GetX();
  double *y = yy->GetY();
  int yn=yy->GetN();
  int xn=xx->GetN();
  int lengthx=xn;
  int lengthy=yn;
  int length=0, d, i, n=0;
  vector<double> out, outx, outy;
  double num, ynum, xdenom, ydenom, denom;
  double timescale = time[1]-time[0];

  length=lengthx<=lengthy?xn:yn;

  max_delay=max_delay>length?length:max_delay;

  double mx=0;
  double my=0;
  int throwaway=.1*max_delay;
	
  for(d=-max_delay;d<max_delay;d++){
    num=0.;
    xdenom=0.;
    ydenom=0.;
    //    for(i=windowlow;i<windowhigh;i++){
    for(i=0;i<length;i++){
      if((i+d)<0 | (i+d)>=length){
	continue;
      }
      else{
	num+=(x[i]-mx)*(y[i+d]-my);
	xdenom+=pow(x[i]-mx, 2);
	ydenom+=pow(y[i+d]-my, 2);
      }
    }
    //square the output. finds strongest correlation/anticorrelation.
    //out.push_back(pow(num/sqrt(xdenom*ydenom), 2));
    out.push_back(num/sqrt(xdenom*ydenom));
    outx.push_back(time[(length/2)+d]);
    n++;
    //out.push_back(pow(num/10000, 2));
		
  }

  double maxIndex=TMath::LocMax(out.size(), &out[0]);
  double offset=(maxIndex-(double)max_delay)*timescale;

  if(xcorr_graph==1){
    TGraph *outt = new TGraph(outx.size(), &outx[0], &out[0]);
    return outt;
  }
  //	if(xcorr_graph==2){
  // return 
  //else
  outx.clear();
  for(int i=0;i<xn;i++){
    outx.push_back(time[i]-offset);
    outy.push_back(y[i]);
  }
  TGraph *outt = new TGraph(outx.size(), &outx[0], &outy[0]);
  return outt;
}

TGraph * crossCorrelateWindowed(TGraph * xx, TGraph * yy, TGraph *wingraph, int max_delay=100000){
  double *window=wingraph->GetY();
  double *x = xx->GetY();
  double *time=xx->GetX();
  double *y = yy->GetY();
  int yn=yy->GetN();
  int xn=xx->GetN();
 
  int lengthx=xn;
  int lengthy=yn;
  int length=0, d, i, n=0;
  vector<double> out, outx, outy;
  double num, ynum, xdenom, ydenom, denom;
  double timescale = time[1]-time[0];

  length=lengthx<=lengthy?xn:yn;
  int throwaway=length/10;//throw away highest delays
  max_delay=max_delay>length?length-throwaway:max_delay;

  double mx=0;
  double my=0;

  for(d=-max_delay;d<max_delay;d++){
    num=0.;
    xdenom=1.;
    ydenom=1.;
    for(i=0;i<length;i++){
      if((i+d)<0 | (i+d)>=length){
	continue;
      }
      else{

	double val=(x[i]-mx)*(y[i+d]-my)*window[i];
	num+=val;
	xdenom+=pow(x[i]-mx, 2)*window[i]; 
	ydenom+=pow(y[i+d]-my, 2)*window[i]; 

      }
	       
    }

    out.push_back(num/sqrt(xdenom*ydenom));
    outx.push_back((d)*timescale);
    n++;
  }

  double maxIndex=TMath::LocMax(out.size(), &out[0]);
  double offset=(maxIndex-(double)max_delay)*timescale;

  TGraph *outt = new TGraph(outx.size(), &outx[0], &out[0]);
  return outt;
}




/*************some plotting things****************/




void TUtil::setWarmPalette(){

  const Int_t rgb = 3;
  const Int_t N = 255;

  Double_t stops[rgb] = {0.34, 0.61, 0.84};
  Double_t red[rgb]   = {0.00, 1., 0.90};
  Double_t green[rgb] = {0., 0.0, 1.0};
  Double_t blue[rgb]  = {0.00, 0.0, 0.00};
  TColor::CreateGradientColorTable(rgb, stops, red, green, blue, N);
  gStyle->SetNumberContours(N);

}

void TUtil::setCoolPalette(){

  const Int_t rgb = 3;
  const Int_t N = 255;


  Double_t red[]    = {0., .0, .0, 1., 1.0};
  Double_t green[]  = {0., .1, .9, .0, 1.0};
  Double_t blue[]   = {0., .80, .90, 0.20, 1.0};
  Double_t stops[] = {0., .25, .50, .75, 1.0};
  TColor::CreateGradientColorTable(rgb, stops, red, green, blue, N);
  gStyle->SetNumberContours(N);

}

