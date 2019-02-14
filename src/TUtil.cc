/*
copyright s. prohira and the T576 collaboration 2018
released under the GNU General Public License version 3
*/

#include "TUtil.hh"

static int fN=1;
static int fNSpec=1.;
static int fNi=1;
static TVirtualFFT * fftr2c=TVirtualFFT::FFT(1, &fN, "R2C ES");
static TVirtualFFT * fftr2cSpec=TVirtualFFT::FFT(1, &fNSpec, "R2C ES");
static TVirtualFFT * fftc2r=TVirtualFFT::FFT(1, &fN, "C2R ES");
static TGraph2D *fXfrmGr2D = new TGraph2D();
static TGraph *fXfrmGr = new TGraph();
static TGraph *fPSDGr = new TGraph();


/***************FFT things************************

the functions here utilise the above static members for the FFT.
FFTW calculates the fastest and most optimal way to do
an FFT the first time it is initialized with new parameters. after that,
each subsequent fft of the same size is very fast. 




**************************************************/

// static void init(){

// }

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

TGraph * TUtil::FFT::psd(TGraph * inGr, double rBW, int dbFlag){
  auto xfrm=fft(inGr);
  int n=xfrm->GetN();
  //  double norm=1./(double)sqrt(n);
  auto xx=xfrm->GetX();
  auto re=xfrm->GetY();
  auto im=xfrm->GetZ();
  double yy[n];

  //resolution bandwidth defaults to nyquist
  rBW=rBW==0?xx[n-1]/2.:rBW;
  
  yy[0]=dbFlag==1?vToDbmHz(rBW,re[0]):re[0];
  for(int i=1;i<(n+1)/2;i++){
    yy[i]=vToDbmHz(rBW, re[i], im[i]);
    yy[i]=dbFlag==1?vToDbmHz(rBW,re[i], im[i]):re[i]*re[i]+im[i]*im[i];
  }
  yy[n/2]=dbFlag==1?vToDbmHz(rBW,re[n/2], im[n/2]):re[n/2]*re[n/2]+im[n/2]*im[n/2];

  TGraph * outGr=new TGraph((n/2), xx, yy);

  *fPSDGr=*outGr;
  //delete outGr;
  return outGr;//fPSDGr;
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
  //  *fXfrmGr=*out;
  //  delete out;
  return out;//fXfrmGr;
}

// vector<double> phasorTransform(vector<double> vec,double ph=1.){
//   int size=vec.size();
//   auto fftd=FFTtools::doFFT(size, &vec[0]);
//   for(int i=0;i<size/2+1;i++){
//     double im=fftd[i].im;
//     fftd[i].im=(2.-ph)*fftd[i].re;
//     fftd[i].re=-ph*im;
//   }
//   double * thing=FFTtools::doInvFFT(size, fftd);
//   auto outvec=toVector(thing, size);
//   delete[] fftd;
//   delete thing;
//   return outvec;
  
TGraph * TUtil::FFT::plotPhase(TGraph *inGr){
  auto ht=hilbertTransform(inGr);
  auto outGr=new TGraph(inGr->GetN());
  for(int i=0;i<inGr->GetN();i++){
    outGr->SetPoint(i, inGr->GetX()[i], atan(ht->GetY()[i]/inGr->GetY()[i]));
  }
  return outGr;
}

TGraph * TUtil::FFT::phasorTransform(TGraph *inGr){
  auto infft=fft(inGr);
  int n=infft->GetN();

  for(int i=0;i<n;i++){
    double im=infft->GetZ()[i];
    infft->GetZ()[i]=0;//infft->GetY()[i];
    infft->GetY()[i]=infft->GetY()[i]+im;
  }
  auto outGr=ifft(infft);
  return outGr;
}

TH2D* TUtil::FFT::spectrogram(TGraph *gr, Int_t binsize , Int_t overlap, Int_t zero_pad_length, int win_type, int dbFlag){
  Int_t size = gr->GetN();
  double dt=(gr->GetX()[1]-gr->GetX()[0]); 
  double samprate=1./dt;
  double xmax = size/samprate;
  double xmin = 0;
  zero_pad_length=zero_pad_length<=binsize?binsize:zero_pad_length;
  Int_t num_zeros=(zero_pad_length-binsize)/2;
  //  Int_t nbins = size/overlap;
  int nbins=size/(binsize-overlap);
  xmax=((double)nbins*((double)binsize-(double)overlap))/samprate;
  char*timebuff;
  //double samplerate = size/(xmax-xmin);
  double bandwidth = 1e9*samprate;
  // if(zero_pad_length!=fNSpec){
  //   fNSpec=zero_pad_length;
  //   fftr2cSpec=TVirtualFFT::FFT(1, &fNSpec, "R2C P K");
  // }
  TGraph * in=new TGraph(zero_pad_length);
  TGraph * outt=new TGraph(zero_pad_length);

  vector<double> sX, sY, sZ;
  Int_t start = 0;
  //  Int_t j=0;
  //cout<<size<<" "<<nbins<<" "<<zero_pad_length<<" "<<binsize<<" "<<overlap<<" "<<num_zeros<<" "<<xmax*samprate<<endl;
  TH2D *spectrogramHist=new TH2D("outhist", "spectrogram", nbins-1, xmin, xmax, (zero_pad_length), 0, samprate);
  
  for(int i=0;i<nbins;i++){
    if(start+binsize-1>0){
      int sampnum=0;
      for(int j=0;j<zero_pad_length;j++){
	if(j<num_zeros||j>=binsize+num_zeros){
	  in->SetPoint(j, gr->GetX()[j], 0.);
	}
	else if(j>=num_zeros&&j<binsize+num_zeros){
	  if(sampnum+start<size){
	    in->SetPoint(j, gr->GetX()[j], gr->GetY()[sampnum+start]*window(sampnum, binsize, win_type));
	  }
	  else{
	    in->SetPoint(j, gr->GetX()[j], 0.);
	  }
	  sampnum++;
	}
      }

	//   for(int j=0;j<num_zeros;j++){
      // 	//if((j+start)>=size)break;
      // 	//    for(int j = 0;j<=zero_pad_length;j++){
      // 	in->SetPoint(j, gr->GetX()[j], 0.);
      // }
      // for(int j=num_zeros;j<binsize+num_zeros;j++){
      // 	//      if((j+start)>=size)break;
      // 	in->SetPoint(j, gr->GetX()[j], gr->GetY()[j+start]*window(j-num_zeros, binsize, win_type));//
      // }
      // for(int j=binsize+num_zeros;j<zero_pad_length;j++){
      // 	//if((j+start)>=size)break;
      // 	in->SetPoint(j, gr->GetX()[j], 0);
      // }
      
      outt=TUtil::FFT::psd(in, samprate/2., dbFlag);
      
      for(int j = 0;j<outt->GetN();j++){
	Double_t z = outt->GetY()[j];
	if(!isfinite(z))z=0.;
	spectrogramHist->SetBinContent(i,j,z);//dbm/hz
	//sX.push_back(i*binsize*dt);
	//sY.push_back(outt->GetX()[j]);
	//sZ.push_back(z);
      }
      //    cout<<sX.size()<<endl;
    }
      start+=(binsize-overlap);
    }


  spectrogramHist->GetYaxis()->SetRangeUser(0, spectrogramHist->GetYaxis()->GetXmax()/2.1);
  spectrogramHist->GetXaxis()->SetTitle("Time (ns)");
  spectrogramHist->GetYaxis()->SetTitle("Frequency (GHz)");
  if(dbFlag){
    spectrogramHist->GetZaxis()->SetTitle("dBm/Hz");
  }
  else{
    spectrogramHist->GetZaxis()->SetTitle("W");
  }
  spectrogramHist->GetZaxis()->SetTitleOffset(1.5);

  outt->Delete();
  in->Delete();

  // spectrogramGr->Delevte();
  return spectrogramHist;
  //return outdat;
}


TH2D * TUtil::FFT::avgSpectrograms(vector<TH2D*> inh){
  TH2D *out = (TH2D*)inh[0]->Clone();
  for(int i=1;i<inh.size();i++){
    if(out->GetNbinsX()==inh[i]->GetNbinsX()&&out->GetNbinsY()==inh[i]->GetNbinsY()){
      out->Add(inh[i]);
    }
  }
  out->Scale(1./inh.size());
  return out;
}



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
  D=D==1?N:D;
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
  return vec*=1./((double)sizey);
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
  //  x=x>V.GetNrows()?V.GetNrows():x;
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

TGraph * TUtil::SVD::expandInBasis(TGraph * G, TMatrixD B, int num){
  TVectorD V = toVector(G);
  double dt = G->GetX()[10]-G->GetX()[9];
  return toGraph(expandInBasis(V, B, num), 1./dt, G->GetX()[0]);
}

TGraph * TUtil::SVD::filter(TGraph *G, TMatrixD B, int num){
  TVectorD V = toVector(G);
  auto filter=expandInBasis(V, B, num);
  //  cout<<V.GetNrows()<<" "<<B.GetNcols()<<endl;
  auto oV=V-filter;
  double dt = G->GetX()[10]-G->GetX()[9];
  return toGraph(oV, 1./dt, G->GetX()[0]);
  
}

TVectorD TUtil::SVD::filter(TVectorD V, TMatrixD B, int num){
  auto filter=expandInBasis(V, B, num);
  return V-filter;
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
  int sizex=M.GetNrows();
  int sizey=M.GetNcols();
  //  cout<<sizey<<" "<<sizex<<endl;
  TH2F * map=new TH2F(name, name, sizex, 0.,sizex, sizey, 0,sizey);
  for(int i=0;i<sizex;i++){
    for(int j=0;j<sizey;j++){
      map->SetBinContent(i,j, M[i][j]);
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


double TUtil::vToDbmHz(double bandwidthGSs, double re, double im){
  double val=re*re+im*im;
  return (10.*log10(val/50.))+30-(10.*log10(bandwidthGSs*1.e9));
}

TGraph * TUtil::normalize(TGraph * inGr){
    double length=inGr->GetN();
    double norm=0;
    for(int i=0;i<length;i++){
      norm+=inGr->GetY()[i]*inGr->GetY()[i];
    }
    TGraph *og = (TGraph*)inGr->Clone();
    for(int i=0;i<length;i++)og->GetY()[i]/=sqrt(norm);
    return og;
  }


TGraph2D * TUtil::normalize(TGraph2D * inGr){
    double length=inGr->GetN();
    double norm=0;
    for(int i=0;i<length;i++){
      norm+=inGr->GetZ()[i]*inGr->GetZ()[i];
    }
    TGraph2D *og = (TGraph2D*)inGr->Clone();
    for(int i=0;i<length;i++)og->GetZ()[i]/=sqrt(norm);
    return og;
}


TGraph * TUtil::CDF(TGraph * inGr, int normed){
  auto outGr=new TGraph(inGr->GetN());
    double val=0;
  for(int i=0;i<inGr->GetN();i++){
    outGr->SetPoint(i, inGr->GetX()[i], val);
    val+=inGr->GetY()[i];
  }
  if(normed==1){
    auto ynormGr=TUtil::scale(outGr, 1./outGr->GetY()[outGr->GetN()-1]);
    auto xnormGr=TUtil::stretch(ynormGr, 1./outGr->GetX()[outGr->GetN()-1]);
    outGr=(TGraph*)xnormGr->Clone();
    delete(ynormGr);
    delete(xnormGr);
  }

  return outGr;
}

double TUtil::mean(TGraph *gr, double t_low, double t_high){
  auto data=gr->GetY();
  int N=gr->GetN();

  t_low=t_low<gr->GetX()[0]?gr->GetX()[0]:t_low;
  t_high=t_high>gr->GetX()[N-1]?gr->GetX()[N-1]:t_high;

  double m=0.;

  double n=0.;
  for(int i=0;i<N;i++){
    if(gr->GetX()[i]>t_low&&gr->GetX()[i]<t_high){
      m+=data[i];
      n+=1.;
    }
  }
  
  return m/n;
}

TGraph * TUtil::removeMean(TGraph *gr, double t_low, double t_high){
  double m=TUtil::mean(gr, t_low, t_high);
  return shiftY(gr, -m);
}

int TUtil::removeMeanInPlace(TGraph *gr, double t_low, double t_high){
  double m=TUtil::mean(gr, t_low, t_high);
  for(int i=0;i<gr->GetN();i++){
    gr->GetY()[i]+=(-m);
  }
  return 1;
}


double TUtil::integrate(TGraph * gr, double t_low, double t_high){
  t_low=t_low>0.?t_low:0.;
  t_high>gr->GetX()[gr->GetN()-1]?gr->GetX()[gr->GetN()-1]:t_high;
  double dt = gr->GetX()[1]-gr->GetX()[0];
  double integral=0.;
  for(int i=0;i<gr->GetN();i++){
    if(gr->GetX()[i]<t_low||gr->GetX()[i]>t_high)continue;
    integral+=gr->GetY()[i]*dt;
  }
  return integral;
}

double TUtil::integratePower(TGraph * gr, double t_low, double t_high){
  t_low=t_low>0.?t_low:0.;
  t_high>gr->GetX()[gr->GetN()-1]?gr->GetX()[gr->GetN()-1]:t_high;
  double dt = gr->GetX()[1]-gr->GetX()[0];
  double integral=0.;
  for(int i=0;i<gr->GetN();i++){
    if(gr->GetX()[i]<t_low||gr->GetX()[i]>t_high)continue;
    integral+=(gr->GetY()[i]*gr->GetY()[i])*dt;
  }
  return integral;
}




TGraph * TUtil::integrateByBin(TGraph *gr, double binNS){
  int nbins=gr->GetX()[gr->GetN()-1]/binNS;
  auto outGr=new TGraph(nbins);

  double t=binNS;
  for(int i=0;i<nbins;i++){
    outGr->SetPoint(i, t, integrate(gr, t-binNS, t));
    t+=binNS;
  }

  return outGr;

}

TGraph * TUtil::avgGraph(vector<TGraph*> inGr){
  auto outGr=new TGraph(inGr[0]->GetN());
  for(int i=0;i<inGr.size();i++){
    for(int j=0;j<outGr->GetN();j++){
      outGr->SetPoint(j, inGr[0]->GetX()[j], inGr[i]->GetY()[j]+outGr->GetY()[j]);
    }
  }
    return TUtil::scale(outGr, 1./inGr.size());
}

TGraph * TUtil::add(TGraph *g1, TGraph *g2, double constant){

  int len=g1->GetN()<g2->GetN()?g1->GetN():g2->GetN();
  TGraph *outGr=new TGraph(len);  
  for(int i=0;i<len;i++){
    outGr->SetPoint(i, g1->GetX()[i], g1->GetY()[i]+(constant*g2->Eval(g1->GetX()[i])));
  }
  outGr->GetXaxis()->SetTitle(g1->GetXaxis()->GetTitle());
  outGr->GetYaxis()->SetTitle(g1->GetYaxis()->GetTitle());

  return outGr;
}

TGraph2D * TUtil::add(TGraph2D *g1, TGraph2D *g2, double constant){

  int len=g1->GetN()<g2->GetN()?g1->GetN():g2->GetN();
  TGraph2D *outGr=new TGraph2D(len);  
  for(int i=0;i<len;i++){
    outGr->SetPoint(i, g1->GetX()[i], g1->GetY()[i], g1->GetZ()[i]+(g2->GetZ()[i]*constant));
  }
  outGr->GetXaxis()->SetTitle(g1->GetXaxis()->GetTitle());
  outGr->GetYaxis()->SetTitle(g1->GetYaxis()->GetTitle());

  return outGr;
}

double TUtil::dot(TGraph *g1, TGraph *g2){
  auto x = g1->GetY();
  auto y = g2->GetY();
  double num=0., xdenom=0., ydenom=0.;
  auto n=g1->GetN()<g2->GetN()?g1->GetN():g2->GetN();
  for(int i=0;i<n;i++){
    num+=(x[i])*(y[i]);
    xdenom+=pow(x[i], 2);
    ydenom+=pow(y[i], 2);
  }
  return num/sqrt(xdenom*ydenom);
  
}

TGraph * TUtil::scale(TGraph *g1, double factor){
  TGraph *outGr=new TGraph(g1->GetN());
  for(int i=0;i<g1->GetN();i++){
    outGr->SetPoint(i, g1->GetX()[i], g1->GetY()[i]*factor);
  }
  outGr->SetName(g1->GetName());
  outGr->SetTitle(g1->GetTitle());
  outGr->GetXaxis()->SetTitle(g1->GetXaxis()->GetTitle());
  outGr->GetYaxis()->SetTitle(g1->GetYaxis()->GetTitle());

  return outGr;
  
}

TGraph * TUtil::stretch(TGraph *g1, double factor){
  TGraph *outGr=new TGraph(g1->GetN());
  for(int i=0;i<g1->GetN();i++){
    outGr->SetPoint(i, g1->GetX()[i]*factor, g1->GetY()[i]);
  }
  outGr->SetName(g1->GetName());
  outGr->SetTitle(g1->GetTitle());
  outGr->GetXaxis()->SetTitle(g1->GetXaxis()->GetTitle());
  outGr->GetYaxis()->SetTitle(g1->GetYaxis()->GetTitle());

  return outGr;

}


TGraph * TUtil::shiftY(TGraph *g1, double factor){
  TGraph *outGr=new TGraph(g1->GetN());
  for(int i=0;i<g1->GetN();i++){
    outGr->SetPoint(i, g1->GetX()[i], g1->GetY()[i]+factor);
  }
  outGr->SetName(g1->GetName());
  outGr->SetTitle(g1->GetTitle());
  outGr->GetXaxis()->SetTitle(g1->GetXaxis()->GetTitle());
  outGr->GetYaxis()->SetTitle(g1->GetYaxis()->GetTitle());

  return outGr;
  
}
TGraph * TUtil::shiftX(TGraph *g1, double factor){
  TGraph *outGr=new TGraph(g1->GetN());
  for(int i=0;i<g1->GetN();i++){
    outGr->SetPoint(i, g1->GetX()[i]+factor, g1->GetY()[i]);
  }
  outGr->SetName(g1->GetName());
  outGr->SetTitle(g1->GetTitle());
  outGr->GetXaxis()->SetTitle(g1->GetXaxis()->GetTitle());
  outGr->GetYaxis()->SetTitle(g1->GetYaxis()->GetTitle());

  return outGr;
  
}




double TUtil::sinc(double x){
  if(x==0.){
    return 1.;
  }
  return sin(pi* x)/(pi*x);
}


TGraph * TUtil::sincInterpolateGraph(TGraph *inGr, double interpGSs){
  double T = inGr->GetX()[1]-inGr->GetX()[0];
  double dt = 1./interpGSs;
  vector<double> xx, yy;
  double t=0;
  double num=0.;
  
  t=0;    
  while(t<inGr->GetX()[inGr->GetN()-1]){
    xx.push_back(t);
    double temp=0;
    for(int i=0;i<inGr->GetN();i++){
      temp+=inGr->GetY()[i]*TUtil::sinc((t-((double)i*T))/T);

    }
  
    yy.push_back(temp);
    t+=dt;
  }
  TGraph *outGr = new TGraph(xx.size(), &xx[0], &yy[0]);
  return outGr;
}


// TGraph * TUtil::sincInterpolateGraphSimple(TGraph *inGr, double interpGSs){
//   double T = inGr->GetX()[1]-inGr->GetX()[0];
//   double dt = 1./interpGSs;
//   vector<double> xx, yy;
//   double t=0;
//   double num=0.;
//   double t_ex = inGr->GetX()[inGr->GetN()/2];
//   double ex=0.;
//   for(int i=0;i<inGr->GetN();i++){
//     ex+=inGr->GetY()[i]*TUtil::sinc((t_ex-((double)i*T))/T);
//   }
  
//   t=0;
//   int i=0;
//   while(t<inGr->GetX()[inGr->GetN()-1]){
//     xx.push_back(t);
//     i=fmod(t, inGr->GetY()[i])>(T/2)?i+1:i;
//     double temp=inGr->GetY()[i]*TUtil::sinc(fmod(t, T)/T);
//     cout<<t<<" "<<i<<" "<<inGr->GetY()[i]<<endl;
//     yy.push_back(temp);
//     t+=dt;
//   }
//   TGraph *outGr = new TGraph(xx.size(), &xx[0], &yy[0]);
//   return outGr;
// }


int TUtil::getInterpolatedGraph(TGraph * inGraph, TGraph *outGraph, double interpGSs){
  ROOT::Math::Interpolator interp(inGraph->GetN(), ROOT::Math::Interpolation::kAKIMA);
  interp.SetData(inGraph->GetN(), inGraph->GetX(), inGraph->GetY());

  //get dt, assuming even sampling.
  double inDt=inGraph->GetX()[50]-inGraph->GetX()[49];
  double inGSs=1./inDt;

  double outDt=1./interpGSs;

  int samps=(int) (inGraph->GetN()*(interpGSs/inGSs));

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

TGraph * TUtil::getChunkOfGraph(TGraph *ingr, double start, double end, int delay_to_zero){
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
  outg->SetTitle(ingr->GetTitle());
  outg->SetName(ingr->GetName());
  outg->GetXaxis()->SetTitle(ingr->GetXaxis()->GetTitle());
  outg->GetYaxis()->SetTitle(ingr->GetYaxis()->GetTitle());
  if(delay_to_zero==0){
    return outg;
  }
  return delayGraph(outg, -start);

}

// TGraph * TUtil::getChunkOfGraph(TGraph *ingr, double start, double end){
//   double *xx=ingr->GetX();
//   double *yy=ingr->GetY();
//   vector<double> outx, outy;
//   double xincr=xx[10]-xx[9];
//   double time=start;
//   while(time<end){
//     outx.push_back(time);
//     outy.push_back(ingr->Eval(time));
//     time+=xincr;
//   }
  
//   TGraph * outg=new TGraph(outx.size(), &outx[0], &outy[0]);
//   return outg;
// }



TGraph * TUtil::delayGraph(TGraph *ingr, double delay){
  double*xx=ingr->GetX();
  double xxout[ingr->GetN()];
  double*yy=ingr->GetY();
  for(int i=0;i<ingr->GetN();i++){
    xxout[i]=(double)xx[i]+delay;
  }
  TGraph *dg=new TGraph(ingr->GetN(), xxout, yy);
  dg->SetTitle(ingr->GetTitle());
  dg->SetName(ingr->GetName());
  dg->GetXaxis()->SetTitle(ingr->GetXaxis()->GetTitle());
  dg->GetYaxis()->SetTitle(ingr->GetYaxis()->GetTitle());

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
  dg->SetTitle(ingr->GetTitle());
  dg->SetName(ingr->GetName());
  dg->GetXaxis()->SetTitle(ingr->GetXaxis()->GetTitle());
  dg->GetYaxis()->SetTitle(ingr->GetYaxis()->GetTitle());

  *outgr=*dg;
  delete dg;
  return 1;
}


TH1F * TUtil::plotResiduals(TGraph *gr1, TGraph *gr2, int nbins, double min, double max){
  TH1F *hist =new TH1F("", "", nbins,min, max);
  for(int i=0;i<gr1->GetN();i++)hist->Fill(gr1->GetY()[i]-gr2->GetY()[i]);
  return hist;
}


TGraph * TUtil::crossCorrelate(TGraph * gr1, TGraph * gr2, double max_delay, double t_low, double t_high){
  double *x = gr1->GetY();
  double *time=gr1->GetX();
  double *y = gr2->GetY();
  int yn=gr1->GetN();
  int xn=gr2->GetN();

  int lengthx=xn;
  int lengthy=yn;
  int length=0;
  vector<double> out, outx, outy;
  double num, ynum, xdenom, ydenom, denom;
  double timescale = time[1]-time[0];

  length=lengthx<=lengthy?xn:yn;
  length=lengthx<=lengthy?xn:yn;
  double throwaway=time[xn-1]/10.;//throw away highest delays, they are unstable
  max_delay=max_delay>time[xn-1]?time[xn-1]-throwaway:max_delay;


  t_high=t_high>=time[xn-1]?time[xn-1]:t_high;
  t_low=t_low<time[0]?time[0]:t_low;

  int max_delay_index=(1./timescale)*max_delay;
  double mx=0;
  double my=0;

  int n=0;
  double t=-max_delay;
  for(int n=-max_delay_index;n<max_delay_index;n++){
    //if(time[d]>=-max_delay&&time[d]<=max_delay){
      num=0.;
      xdenom=0.;
      ydenom=0.;
      for(int i=0;i<length;i++){
	if((i+n)>0 && (i+n)<length && time[i]>=t_low && time[i]<=t_high){
	  
	  num+=(x[i]-mx)*(y[i+n]-my);
	  xdenom+=pow(x[i]-mx, 2);
	  ydenom+=pow(y[i+n]-my, 2);
	  
	}
      }
      out.push_back(num/sqrt(xdenom*ydenom));
      //      outx.push_back(time[(length/2)+n]);
      outx.push_back((double)n *timescale);
      //    n++;    
    }



  
  TGraph *outt = new TGraph(outx.size(), &outx[0], &out[0]);
  outt->GetXaxis()->SetTitle("offset (ns)");
  outt->GetYaxis()->SetTitle("CC coefficient");
  outt->GetYaxis()->SetTitleOffset(1.15);
  outt->SetTitle("");
  outt->GetXaxis()->SetRangeUser(-max_delay, max_delay);
  return outt;
}


TGraph * TUtil::crossCorrelateWindowed(TGraph * gr1, TGraph * gr2, TGraph *grWindow, double max_delay, double t_low, double t_high){
  double *x = gr1->GetY();
  double *time=gr1->GetX();
  double *y = gr2->GetY();
  double *window = grWindow->GetY();
  int yn=gr1->GetN();
  int xn=gr2->GetN();

  int lengthx=xn;
  int lengthy=yn;
  int length=0;
  vector<double> out, outx, outy;
  double num, ynum, xdenom, ydenom, denom;
  double timescale = time[1]-time[0];

  length=lengthx<=lengthy?xn:yn;
  length=lengthx<=lengthy?xn:yn;
  double throwaway=time[xn-1]/10.;//throw away highest delays, they are unstable
  max_delay=max_delay>time[xn-1]?time[xn-1]-throwaway:max_delay;

  t_high=t_high>=time[xn-1]?time[xn-1]:t_high;
  t_low=t_low<0.?0.:t_low;

  int max_delay_index=(1./timescale)*max_delay;
  double mx=0;
  double my=0;

  int n=0;
  double t=-max_delay;
  for(int n=-max_delay_index;n<max_delay_index;n++){
    //if(time[d]>=-max_delay&&time[d]<=max_delay){
      num=0.;
      xdenom=0.;
      ydenom=0.;
      for(int i=0;i<length;i++){
	if((i+n)>0 && (i+n)<length && time[i]>=t_low && time[i]<=t_high){
	  
	  num+=(x[i]-mx)*(y[i+n]-my);
	  xdenom+=pow(x[i]-mx, 2);
	  ydenom+=pow(y[i+n]-my, 2);
	  
	}
      }
      out.push_back(num/sqrt(xdenom*ydenom));
      outx.push_back(time[(length/2)+n]);
      //    n++;    
    }


  
  TGraph *outt = new TGraph(outx.size(), &outx[0], &out[0]);

  return outt;
}


TGraph * TUtil::align(TGraph * gr1, TGraph * gr2, double max_delay, double t_low, double t_high){
  double *x = gr1->GetY();
  double *time=gr1->GetX();
  double *y = gr2->GetY();
  int yn=gr1->GetN();
  int xn=gr2->GetN();
  int lengthx=xn;
  int lengthy=yn;
  int length=0;
  vector<double> out, outx, outy;
  double num, ynum, xdenom, ydenom, denom;
  double timescale = time[1]-time[0];

  length=lengthx<=lengthy?xn:yn;
  length=lengthx<=lengthy?xn:yn;

  double throwaway=time[xn-1]/10.;//throw away highest delays, they are unstable
  max_delay=max_delay>time[xn-1]?time[xn-1]-throwaway:max_delay;

  t_high=t_high>=time[xn-1]?time[xn-1]:t_high;
  t_low=t_low<time[0]?time[0]:t_low;

  int max_delay_index=(1./timescale)*max_delay;
  double mx=0;
  double my=0;

  int n=0;
  double t=-max_delay;
  for(int n=-max_delay_index;n<max_delay_index;n++){
    //if(time[d]>=-max_delay&&time[d]<=max_delay){
      num=0.;
      xdenom=0.;
      ydenom=0.;
      for(int i=0;i<length;i++){
	if((i+n)>0 && (i+n)<length && time[i]>=t_low && time[i]<=t_high){
	  
	  num+=(x[i]-mx)*(y[i+n]-my);
	  xdenom+=pow(x[i]-mx, 2);
	  ydenom+=pow(y[i+n]-my, 2);
	  
	}
      }
      out.push_back(num/sqrt(xdenom*ydenom));
      outx.push_back(time[(length/2)+n]);
      //    n++;    
    }

  double maxIndex=TMath::LocMax(out.size(), &out[0]);
  double offset=(maxIndex-(double)max_delay_index)*timescale;

  // outx.clear();
  // for(int i=0;i<xn;i++){
  //   outx.push_back(time[i]-offset);
  //   outy.push_back(y[i]);
  // }

  // TGraph *outt = new TGraph(outx.size(), &outx[0], &outy[0]);
  
  return delayGraph(gr2, -offset);
}


TGraph * TUtil::alignToOther(TGraph * gr1, TGraph * gr2, TGraph *othGr, double max_delay, double t_low, double t_high){
  double *x = gr1->GetY();
  double *time=gr1->GetX();
  double *y = gr2->GetY();
  double *yOth=othGr->GetY();
  int yn=gr1->GetN();
  int xn=gr2->GetN();
  int lengthx=xn;
  int lengthy=yn;
  int length=0;
  vector<double> out, outx, outy;
  double num, ynum, xdenom, ydenom, denom;
  double timescale = time[1]-time[0];

  length=lengthx<=lengthy?xn:yn;
  length=lengthx<=lengthy?xn:yn;

  double throwaway=time[xn-1]/10.;//throw away highest delays, they are unstable
  max_delay=max_delay>time[xn-1]?time[xn-1]-throwaway:max_delay;

  t_high=t_high>=time[xn-1]?time[xn-1]:t_high;
  t_low=t_low<time[0]?time[0]:t_low;

  int max_delay_index=(1./timescale)*max_delay;
  double mx=0;
  double my=0;

  int n=0;
  double t=-max_delay;
  for(int n=-max_delay_index;n<max_delay_index;n++){
    //if(time[d]>=-max_delay&&time[d]<=max_delay){
      num=0.;
      xdenom=0.;
      ydenom=0.;
      for(int i=0;i<length;i++){
	if((i+n)>0 && (i+n)<length && time[i]>=t_low && time[i]<=t_high){
	  
	  num+=(x[i]-mx)*(y[i+n]-my);
	  xdenom+=pow(x[i]-mx, 2);
	  ydenom+=pow(y[i+n]-my, 2);
	  
	}
      }
      out.push_back(num/sqrt(xdenom*ydenom));
      outx.push_back(time[(length/2)+n]);
      //    n++;    
    }

  double maxIndex=TMath::LocMax(out.size(), &out[0]);
  double offset=(maxIndex-(double)max_delay_index)*timescale;

  outx.clear();
  for(int i=0;i<othGr->GetN();i++){
    outx.push_back(othGr->GetX()[i]-offset);
    //    outy.push_back(y[i]);
  }

  TGraph *outt = new TGraph(outx.size(), &outx[0], yOth);
  
  return outt;
}


vector<TGraph*> TUtil::alignMultiple(vector<TGraph*> inGr, double max_delay, double t_low, double t_high){   
  vector<TGraph*>outgraphs;
  TGraph *g1=inGr[0];
  //  g1->Draw("al PLC");
  outgraphs.push_back(g1);
  for(int i=1;i<inGr.size();i++){
    outgraphs.push_back(align(g1, inGr[i], max_delay, t_low, t_high));
    //    cout<<i<<endl;
    //    outgraphs[i]->Draw("l same PLC");
  }
  return outgraphs;
}

vector<TGraph*> TUtil::alignMultipleToOther(vector<TGraph*> inGr, vector<TGraph*> othGr, double max_delay, double t_low, double t_high){   
  vector<TGraph*>outgraphs;
  TGraph *g1=inGr[0];
  //  g1->Draw("al PLC");
  outgraphs.push_back(othGr[0]);
  for(int i=1;i<inGr.size();i++){
    outgraphs.push_back(alignToOther(g1, inGr[i], othGr[i], max_delay, t_low, t_high));
    //    cout<<i<<endl;
    //    outgraphs[i]->Draw("l same PLC");
  }
  return outgraphs;
}


TGraph * TUtil::makeCW(double freq,  double amp, double t_min, double t_max, double GSs, double phase){

  int n=(t_max-t_min)*GSs;
  TGraph * oG=new TGraph(n);
  double dt=1./GSs;
  for(int i=0;i<n;i++){
    double t = (double)i*dt;
    double temp=amp*sin(2.*pi*freq*t + phase);
    oG->SetPoint(i, t, temp);
  }
  return oG;
}


TGraph * TUtil::lowpassFilter(TGraph *ingr, double cutoff, int order){
  double * yy=ingr->GetY();
  double *xx= ingr->GetX();
  int n = ingr->GetN();
  vector<double> outx, outy;
  
  double w = cutoff*2.*pi*1.e9;
  double T = (xx[10]-xx[9])*1.e-9;;
  double a, b, c, value;
  
  //cout<<setprecision(12);
  //cout<<"filter coefficients"<<endl<<a<<endl<<b<<endl<<c<<endl;
  //  int size = in.size();
  if(order==1){
    //    a = w*T;
    b = exp(-w*T);
    a=1.-b;
    
    for(int i=0;i<n;i++){
      if(i>0){
	
	value = a*yy[i]+b*outy[i-1];
	
	outy.push_back(value);
      }
      if(i==0){
	value = a*yy[i];
	//value=0.;
	outy.push_back(value);
      } 
    }
  }
  
  else{
    //a = pow(T, 2.)*pow(w, 2.)*exp(-2.*w*T);
    b = 2.*exp(-w*T);
    c = exp(-2.*w*T);
    a=1.-(b-c);	
    //	cout<<"filter coefficients"<<endl<<a<<endl<<b<<endl<<c<<endl;

    for(int i=0;i<n;i++){
      if(i>1){
	value = a*yy[i-1]+b*outy[i-1]-c*outy[i-2];
	outy.push_back(value);

      }
      if(i==1){
	value = a*yy[i-1]+b*outy[i-1];
	outy.push_back(value);
			
      }
      if(i==0){
	outy.push_back(0.);
      } 
      // outx.push_back(xx[i]);
    }
  }
  TGraph *outgr = new TGraph(n, ingr->GetX(), &outy[0]);
  return outgr;
}

TGraph * TUtil::brickWallFilter(TGraph * inGr, double low, double high){
  
  auto fT=TUtil::FFT::fft(inGr);

  auto fs=inGr->GetX()[1]-inGr->GetX()[0];
  auto df=fs/inGr->GetN();

  int indL=low/df;
  int indH=high/df;

  for(int i=0;i<indL;i++){
    fT->SetPoint(i, fT->GetX()[i], 0., 0.);
  }

  for(int i=indH;i<inGr->GetN();i++){
    fT->SetPoint(i, fT->GetX()[i], 0., 0.);
  }

  auto outGr=(TGraph*)TUtil::FFT::ifft(fT)->Clone();
  return outGr;
}



TGraph * TUtil::makeNullData(TGraph *sig, TGraph *back, double t_min, double t_max, double scale){
  auto sigchunk=getChunkOfGraph(sig, 0., (t_max-t_min));
  auto backchunk=getChunkOfGraph(back, t_min, t_max ,1);
  return add(sigchunk, TUtil::scale(backchunk, scale));
}


double TUtil::integrate2D(TH2D *h, double xmin, double xmax, double ymin, double ymax, double & err){
  Int_t xmin_bin = h->GetXaxis()->FindBin(xmin);
  Int_t xmax_bin = h->GetXaxis()->FindBin(xmax);
  Int_t ymin_bin = h->GetYaxis()->FindBin(ymin);
  Int_t ymax_bin = h->GetYaxis()->FindBin(ymax);

  return  h->IntegralAndError(xmin_bin, xmax_bin, ymin_bin, ymax_bin, err);
  
}


// double integral_1d(TH1F *h, double xmin, double xmax, double & err){
//   Int_t xmin_bin = h->GetXaxis()->FindBin(xmin);
//   Int_t xmax_bin = h->GetXaxis()->FindBin(xmax);
 

//   return  h->IntegralAndError(xmin_bin, xmax_bin, err);
  
// }

/* 
here is the layout for the coordinates used. centermost band is the 
signal band. those to either side in both dims are the sidebands. the integrals
ix1, iy2 etc are the integrals of those quadrants. ib1-4 are averaged
for the overall background. ix1-background and ix2-background are averaged
to get the signal band background, same in y. finally, background and both 
signal band backgrounds are subtracted from the signal quadrant to get signal. 
   __________________________
  |      |   |   |   |       |
  |______|___|___|___|_______|y4 
  |      |ib3|iy2|ib4|       | 
  |______|___|___|___|_______|y3
  |      |ix1|sig|ix2|       | 
  |______|___|___|___|_______|y2
  |      |ib1|iy1|ib2|       | 
  |______|___|___|___|_______|y1
  |      |   |   |   |       | 
  |      |   |   |   |       | 
  |______|___|___|___|_______| 
        x1   x2  x3  x4 
    


 */


double TUtil::sidebandSubtraction2D(TH2D *h, double sband_x1, double sband_x2, double sband_y1, double sband_y2, double & err, int draw){

  double x1, x2, x3, x4, y1, y2, y3, y4, bandwidth_x, bandwidth_y,ix1, ix2, iy1, iy2, isig, ib1, ib2, ib3, ib4, background, bandy, bandx, avg_x, avg_y, sig;
  double ix1_err, ix2_err, iy1_err, iy2_err, isig_err, ib1_err, ib2_err, ib3_err, ib4_err, background_err, bandy_err, bandx_err, avg_x_err, avg_y_err, sig_err;
  //assign the cordinates of the signal and sidebands
  x2 = sband_x1;
  x3 = sband_x2;
  y2 = sband_y1;
  y3 = sband_y2;
  bandwidth_x = x3-x2;
  bandwidth_y = y3-y2;
  x1 = x2-bandwidth_x;
  x4 = x3+bandwidth_x;
  y1 = y2-bandwidth_y;
  y4 = y3+bandwidth_y;
  //make the background integrals and average
  ib1 = integrate2D(h, x1, x2, y1, y2, ib1_err);
  ib2 = integrate2D(h, x3, x4, y1, y2, ib2_err);
  ib3 = integrate2D(h, x1, x2, y3, y4, ib3_err);
  ib4 = integrate2D(h, x3, x4, y3, y4, ib4_err);

  background =(ib1+ib2+ib3+ib4)/4;

  //make the signal band backgrounds
  ix1 = integrate2D(h, x1, x2, y2, y3, ix1_err);
  ix2 = integrate2D(h, x3, x4, y2, y3, ix2_err);
  iy1 = integrate2D(h, x2, x3, y1, y2, iy1_err);
  iy2 = integrate2D(h, x2, x3, y3, y4, iy2_err);

   bandx = ((ix1-background)+(ix2-background))/2.;
   bandy = ((iy1-background)+(iy2-background))/2.;
 
 
  //make the signal integral
  isig = integrate2D(h, x2, x3, y2, y3, isig_err);
  

  //errors. add in quadrature and average.
  background_err = sqrt(pow(ib1_err, 2)+pow(ib2_err, 2)+pow(ib3_err, 2)+pow(ib4_err, 2))/2;
  bandx_err = sqrt(pow(ix1_err, 2)+pow(ix2_err, 2))/sqrt(2);
  bandy_err = sqrt(pow(iy1_err, 2)+pow(iy2_err, 2))/sqrt(2);

  sig_err = sqrt(pow(background_err, 2)+pow(bandx_err, 2)+pow(bandy_err, 2));
  err=sig_err;
  //construct signal
  sig = isig-(background+bandx+bandy);
  //  sig = isig-((ix1+ix2+iy1+iy2)/4);

  //  cout<<"total events: "<<sig<<"+/-"<<sig_err<<endl;
  //  TString tit = std::to_string(sig);

  if(draw==1){
    h->Draw("colz");
    double ymin = gPad->GetUymin();//h->GetYaxis()->GetXmin();
    double ymax = gPad->GetUymax();//h->GetYaxis()->GetXmax();
    double xmin = gPad->GetUxmin();//h->GetXaxis()->GetXmin();
    double xmax = gPad->GetUxmax();//h->GetXaxis()->GetXmax();

    TLine *l1 = new TLine(x2, ymin, x2, ymax);
    TLine *l2 = new TLine(x3, ymin, x3, ymax);
    TLine *l3 = new TLine(xmin, y2, xmax, y2);
    TLine *l4 = new TLine(xmin, y3, xmax, y3);
    TLine *l5 = new TLine(x1, ymin, x1, ymax);
    TLine *l6 = new TLine(x4, ymin, x4, ymax);
    TLine *l7 = new TLine(xmin, y1, xmax, y1);
    TLine *l8 = new TLine(xmin, y4, xmax, y4);
    l1->SetLineColor(kRed);
    l2->SetLineColor(kRed);
    l3->SetLineColor(kRed);
    l4->SetLineColor(kRed);

    l1->Draw();
    l2->Draw();
    l3->Draw();
    l4->Draw();
    l5->Draw();
    l6->Draw();
    l7->Draw();
    l8->Draw();
    char title[100];

    h->SetStats(0);
    ((TCanvas*)gROOT->GetListOfCanvases()->At(0))->SetRightMargin(0.15);
  }


  
  return sig;
} 

/*
 here is the layout for the coordinates used. centermost band is the 
 signal band. those to either side in both dims are the sidebands. the integrals
 ix1, iy2 etc are the integrals of those quadrants. ib1-4 are averaged
 for the overall background. ix1-background and ix2-background are averaged
 to get the signal band background, same in y. finally, background and both 
 signal band backgrounds are subtracted from the signal quadrant to get signal. 
    __________________________
   |      |   |   |   |       |
   |______|___|___|___|_______|y4 
   |      |b13|y34|b33|       | 
   |______|___|___|___|_______|y3
   |      |x12|sig|x34|       | 
   |______|___|___|___|_______|y2
   |      |b11|y12|b31|       | 
   |______|___|___|___|_______|y1
   |      |   |   |   |       | 
   |      |   |   |   |       | 
   |______|___|___|___|_______| 
         x1   x2  x3  x4 
    
*/

// double sidebandSubtract(TH2D* M, x2, y2, nbins):
//     x1=x2-nbins;
//     y1=y2-nbins;
//     x3=x2+nbins;
//     y3=y2+nbins;
//     x4=x3+nbins;
//     y4=y3+nbins;

//     x12=0.
//     x34=0.
//     y12=0.
//     y34=0.
//     b13=0.
//     b31=0.
//     b11=0.
//     b33=0.
    
//     sig=0.
//     #x12
//     for i in range(x1,x2):
//         for j in range(y2,y3):
//             x12+=M[i][j]
//     #x34
//     for i in range(x3,x4):
//         for j in range(y2,y3):
//             x34+=M[i][j]

//     #y12
//     for i in range(x2,x3):
//         for j in range(y1,y2):
//             y12+=M[i][j]

//     #y34
//     for i in range(x2,x3):
//         for j in range(y3,y4):
//             y34+=M[i][j]

//     #b11
//     for i in range(x1,x2):
//         for j in range(y1,y2):
//             b11+=M[i][j]
//     #b31
//     for i in range(x3,x4):
//         for j in range(y1,y2):
//             b31+=M[i][j]

//     #b13
//     for i in range(x1,x2):
//         for j in range(y3,y4):
//             b13+=M[i][j]

//     #b33
//     for i in range(x3,x4):
//         for j in range(y3,y4):
//             b33+=M[i][j]

    

            
//     #sig
//     for i in range(x2,x3):
//         for j in range(y2,y3):
//             sig+=M[i][j]

//     avgy=(y12+y34)/2.
//     avgx=(x12+x34)/2.
//     bkgnd=(b11+b13+b31+b33)/4.
    
//     #testing
//     #outsig=sig-avgy-avgx+bkgnd;
//     outsig=sig-avgy;

//     return outsig


//these are all from wikipedia.
double TUtil::window(int i, int n, int type){
  switch (type){
    case 0:
      return bartlettWindow(i, n);
    case 1:
      return welchWindow(i, n);
    case 2:
      return hannWindow(i, n);
    case 3:
      return blackmanNuttallWindow(i, n);
    default:
      return 1.;
      break;
    }
}

double TUtil::bartlettWindow(int i, int n){
  return 1.-abs((2.*(double)i-(double)n)/((double)n));
}

double TUtil::welchWindow(int i, int n){
  return 1.-pow((2.*(double)i-(double)n)/((double)n), 2);
}

double TUtil::hannWindow(int i, int n){
  return sin(pi*i/(n))*sin(pi*i/(n));
}
//i have no idea why there are increments in the pi coeff...
double TUtil::blackmanNuttallWindow(int i, int n){
  return .3635819-.4891775*cos(2.*pi*i/(n)) + .1365995*cos(4.*pi*i/(n)) - .0106411*cos(6.*pi*i/(n));
}

double TUtil::deg2Rad(double deg) {
  return (deg * pi / 180.);
}

double TUtil::rad2Deg(double rad) {
  return (rad * 180. / pi);
}

/*************some plotting things****************/


void TUtil::draw(vector<TGraph*> inGr){
  inGr[0]->Draw("al PLC");
  for(int i=0;i<inGr.size();i++){
    inGr[i]->Draw("l same PLC");
  }
}

// void TUtil::draw(TGraph** inGr){
//   inGr[0]->Draw("al PLC");
//   for(int i=0;i<inGr.size();i++){
//     inGr[i]->Draw("l same PLC");
//   }
// }


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

