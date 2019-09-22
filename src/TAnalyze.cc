/*
  copyright s. prohira and the T576 collaboration 2018
  released under the GNU General Public License version 3
*/
#include <t576/TAnalyze.hh>

TH2D* TAnalyze::avgSpectrogram(int ch, int dataset, int major, int minor, int nfft, int nOverlap, int padTo, int window, int log, int norm){
  auto ev=new T576Event(50, dataset);
  auto hist=ev->drawAvgSpectrogram(major, minor, 0,ch, ev->scopeNEvents, nfft, nOverlap, padTo, window, log);
  if(norm==1){
    TUtil::normalize(hist, 0, ev->analogBandwidth);
  }
  delete ev;
  return hist;
}

TH2D* TAnalyze::avgSpectrogram(int ch, int dataset, int major, int minor, double tLow, double tHigh, int nfft, int nOverlap, int padTo, int window, int log, int norm){
  auto ev=new T576Event(50, dataset);
  auto hist=ev->drawAvgSpectrogram(major, minor, 0,ch, ev->scopeNEvents, tLow, tHigh, nfft, nOverlap, padTo, window, log);
  if(norm==1){
    TUtil::normalize(hist, 0, ev->analogBandwidth);
  }
  delete ev;
  return hist;
}

int TAnalyze::drawAvgRealNullWithGeom(int ch, int major, int minor, int nfft, int nOverlap, int padTo, int window, int log, int norm){
    auto can=new TCanvas("", "", 1200, 400);
    can->Divide(3, 0);
    can->cd(1)->SetRightMargin(.15);

    auto ev=new T576Event(50, 1);
    auto evNull=new T576Event(50, 3);
    //    TUtil::setColdPalette();
    ev->scope->ch[0]->Draw("al PLC");

    auto avgSpec=ev->drawAvgSpectrogram(major, minor, 0,ch,ev->scopeNEvents, nfft, nOverlap, padTo, window, log);
    avgSpec->SetTitle("Data CH"+TString::Itoa(ch, 10));

    TUtil::ranges(avgSpec, 0, 90, 0, 3);

    auto zmax=avgSpec->GetMaximum();
    auto zmin=avgSpec->GetMinimum();
    can->Draw();
    auto cursors=TUtil::drawPeakCursorXY(avgSpec, kRed);
    auto minx=38;
    auto maxx=50.;
    auto miny=1.9;
    auto maxy=2.25;
    can->cd(2)->SetRightMargin(.15);

    ev->drawGeom(ch);
    auto geom=(TGraph*)gPad->GetPrimitive("geom");
    geom->SetTitle(Form("%.1f GHz", ev->frequency));
    can->cd(3)->SetRightMargin(.15);



    auto avgSpecN=evNull->drawAvgSpectrogram(major, minor, 0,ch,evNull->scopeNEvents, nfft, nOverlap, padTo, window, log);
    if(norm==1){
     
      auto normD=TUtil::norm(avgSpec, 0, ev->analogBandwidth);
      auto normN=TUtil::norm(avgSpecN, 0, ev->analogBandwidth);
      avgSpecN->Scale(normD/normN);
    }

    avgSpecN->SetTitle("Null CH"+TString::Itoa(ch, 10));
    avgSpecN->GetZaxis()->SetRangeUser(zmin, zmax);
    TUtil::ranges(avgSpecN, 0, 90, 0, 3);
    cursors[0]->Draw();
    cursors[1]->Draw();
    can->Draw();
    return 1;
}


int TAnalyze::drawAvgRealNullFull(int ch, int major, int minor, int nfft, int nOverlap, int padTo, int window, int log, int norm){
    auto can=new TCanvas("", "", 1200, 400);
    can->Divide(3, 0);
    can->cd(1)->SetRightMargin(.15);

    auto ev=new T576Event(50, 1);
    auto evNull=new T576Event(50, 3);
    auto evFull=new T576Event(50, 2);
    //    TUtil::setColdPalette();
    ev->scope->ch[0]->Draw("al PLC");

    auto avgSpec=ev->drawAvgSpectrogram(major, minor, 0,ch,ev->scopeNEvents, nfft, nOverlap, padTo, window, log);
    avgSpec->SetTitle("Data CH"+TString::Itoa(ch, 10));

    TUtil::ranges(avgSpec, 0, 90, 0, 3);

    auto zmax=avgSpec->GetMaximum();
    auto zmin=avgSpec->GetMinimum();
    can->Draw();
    auto cursors=TUtil::drawPeakCursorXY(avgSpec, kRed);
    auto minx=38;
    auto maxx=50.;
    auto miny=1.9;
    auto maxy=2.25;
    can->cd(2)->SetRightMargin(.15);
    auto avgSpecFull=evFull->drawAvgSpectrogram(major, minor, 0,ch,evFull->scopeNEvents, nfft, nOverlap, padTo, window, log);
    avgSpecFull->SetTitle("Full CH"+TString::Itoa(ch, 10));
    //avgSpecN->GetZaxis()->SetRangeUser(zmin, zmax);
    TUtil::ranges(avgSpecFull, 0, 90, 0, 3);
    cursors[0]->Draw();
    cursors[1]->Draw();
    // ev->drawGeom(ch);
    // auto geom=(TGraph*)gPad->GetPrimitive("geom");
    // geom->SetTitle(Form("%.1f GHz", ev->frequency));

    can->cd(3)->SetRightMargin(.15);



    auto avgSpecN=evNull->drawAvgSpectrogram(major, minor, 0,ch,evNull->scopeNEvents, nfft, nOverlap, padTo, window, log);
    if(norm==1){
     
      auto normD=TUtil::norm(avgSpec, 0, ev->analogBandwidth);
      auto normN=TUtil::norm(avgSpecN, 0, ev->analogBandwidth);
      avgSpecN->Scale(normD/normN);
    }

    avgSpecN->SetTitle("Null CH"+TString::Itoa(ch, 10));
    avgSpecN->GetZaxis()->SetRangeUser(zmin, zmax);
    TUtil::ranges(avgSpecN, 0, 90, 0, 3);
    cursors[0]->Draw();
    cursors[1]->Draw();
    can->Draw();
    return 1;
}


TCanvas * TAnalyze::avgRealNull(int ch, int major, int minor, int nfft, int nOverlap, int padTo, int window, int dbFlag, int norm, int cursorType){
    auto can=new TCanvas("", "", 1300, 600);
    can->Divide(2, 0);
    can->cd(1)->SetRightMargin(.15);
    can->cd(1)->SetLeftMargin(.12);
 
    auto ev=new T576Event(50, 1);
    auto evNull=new T576Event(50, 3);
    //TUtil::setColdPalette();
    ev->scope->ch[0]->Draw("al PLC");

    auto avgSpec=ev->drawAvgSpectrogram(major, minor, 0,ch,ev->scopeNEvents, nfft, nOverlap, padTo, window, dbFlag);
    avgSpec->SetTitle("Data CH"+TString::Itoa(ch, 10));
    avgSpec->SetName("data");
    TUtil::ranges(avgSpec, 0, 90, 0, 3);

    auto zmax=avgSpec->GetMaximum();
    auto zmin=avgSpec->GetMinimum();
    can->Draw();
    auto cursors=vector<TLine*>();
    if(cursorType==1){
      cursors=TUtil::drawPeakCursorXY(avgSpec, kRed);
    }
    auto minx=38;
    auto maxx=50.;
    auto miny=1.9;
    auto maxy=2.25;

    can->cd(2)->SetRightMargin(.15);
    can->cd(2)->SetLeftMargin(.12);
    auto avgSpecN=evNull->drawAvgSpectrogram(major, minor, 0,ch,evNull->scopeNEvents, nfft, nOverlap, padTo, window, dbFlag);
    avgSpecN->SetName("null");
    if(norm==1){
     
      auto normD=TUtil::norm(avgSpec, 0, ev->analogBandwidth);
      auto normN=TUtil::norm(avgSpecN, 0, ev->analogBandwidth);
      avgSpecN->Scale(normD/normN);
    }
    avgSpecN->SetTitle("Null CH"+TString::Itoa(ch, 10));
    avgSpecN->GetZaxis()->SetRangeUser(zmin, zmax);
    TUtil::ranges(avgSpecN, 0, 90, 0, 3);
    if(cursorType==1){
      cursors[0]->Draw();
      cursors[1]->Draw();
    }
    can->Draw();
    return can;
}

TCanvas * TAnalyze::avgRealFull(int ch, int major, int minor, int nfft, int nOverlap, int padTo, int window, int log, int cursorType){
    auto can=new TCanvas("", "", 1200, 600);
    can->Divide(2, 0);
    can->cd(1)->SetRightMargin(.15);
    can->cd(1)->SetLeftMargin(.12);
    auto ev=new T576Event(50, 1);
    auto evNull=new T576Event(50, 2);
    //TUtil::setColdPalette();
    ev->scope->ch[0]->Draw("al PLC");

    auto avgSpec=ev->drawAvgSpectrogram(major, minor, 0,ch,ev->scopeNEvents, nfft, nOverlap, padTo, window, log);
    avgSpec->SetTitle("Filtered CH"+TString::Itoa(ch, 10));

    TUtil::ranges(avgSpec, 0, 90, 0, 3);

    auto zmax=avgSpec->GetMaximum();
    auto zmin=avgSpec->GetMinimum();
    can->Draw();
    auto cursors=vector<TLine*>();
    if(cursorType==1){
      cursors=TUtil::drawPeakCursorXY(avgSpec, kRed);
    }
    auto minx=38;
    auto maxx=50.;
    auto miny=1.9;
    auto maxy=2.25;

    can->cd(2)->SetRightMargin(.15);
    can->cd(2)->SetLeftMargin(.12);
    auto avgSpecN=evNull->drawAvgSpectrogram(major, minor, 0,ch,evNull->scopeNEvents, nfft, nOverlap, padTo, window, log);
    avgSpecN->SetTitle("Full CH"+TString::Itoa(ch, 10));
    //    avgSpecN->GetZaxis()->SetRangeUser(zmin, zmax);
    TUtil::ranges(avgSpecN, 0, 90, 0, 3);
    if(cursorType==1){
      cursors[0]->Draw();
      cursors[1]->Draw();
    }
    
    
    return can;
}



int TAnalyze::drawAvgRealFullWithGeom(int ch, int major, int minor, int nfft, int nOverlap, int padTo, int window, int log){
    auto can=new TCanvas("", "", 1200, 400);
    can->Divide(3, 0);
    can->cd(1)->SetRightMargin(.15);

    auto ev=new T576Event(50, 1);
    auto evNull=new T576Event(50, 2);
    //    TUtil::setColdPalette();
    ev->scope->ch[0]->Draw("al PLC");

    auto avgSpec=ev->drawAvgSpectrogram(major, minor, 0,ch,ev->scopeNEvents, nfft, nOverlap, padTo, window, log);
    avgSpec->SetTitle("Filtered CH"+TString::Itoa(ch, 10));

    TUtil::ranges(avgSpec, 0, 90, 0, 3);

    auto zmax=avgSpec->GetMaximum();
    auto zmin=avgSpec->GetMinimum();
    can->Draw();
    auto cursors=TUtil::drawPeakCursorXY(avgSpec, kRed);
    auto minx=38;
    auto maxx=50.;
    auto miny=1.9;
    auto maxy=2.25;
    can->cd(2)->SetRightMargin(.15);

    ev->drawGeom(ch);
    auto geom=(TGraph*)gPad->GetPrimitive("geom");
    geom->SetTitle(Form("%.1f GHz", ev->frequency));
    can->cd(3)->SetRightMargin(.15);
    auto avgSpecN=evNull->drawAvgSpectrogram(major, minor, 0,ch,evNull->scopeNEvents, nfft, nOverlap, padTo, window, log);
    avgSpecN->SetTitle("Full CH"+TString::Itoa(ch, 10));
    //avgSpecN->GetZaxis()->SetRangeUser(zmin, zmax);
    TUtil::ranges(avgSpecN, 0, 90, 0, 3);
    cursors[0]->Draw();
    cursors[1]->Draw();
    can->Draw();
    return 1;
}


int TAnalyze::drawAvgRealFull(int ch, int major, int minor, int nfft, int nOverlap, int padTo, int window, int log){
    auto can=new TCanvas("", "", 1200, 600);
    can->Divide(2, 0);
    can->cd(1)->SetRightMargin(.15);
    can->cd(1)->SetLeftMargin(.12);
    auto ev=new T576Event(50, 1);
    auto evNull=new T576Event(50, 2);
    //TUtil::setColdPalette();
    ev->scope->ch[0]->Draw("al PLC");

    auto avgSpec=ev->drawAvgSpectrogram(major, minor, 0,ch,ev->scopeNEvents, nfft, nOverlap, padTo, window, log);
    avgSpec->SetTitle("Filtered CH"+TString::Itoa(ch, 10));

    TUtil::ranges(avgSpec, 0, 90, 0, 3);

    auto zmax=avgSpec->GetMaximum();
    auto zmin=avgSpec->GetMinimum();
    can->Draw();
    auto cursors=TUtil::drawPeakCursorXY(avgSpec, kRed);
    auto minx=38;
    auto maxx=50.;
    auto miny=1.9;
    auto maxy=2.25;

    can->cd(2)->SetRightMargin(.15);
    can->cd(2)->SetLeftMargin(.12);
    auto avgSpecN=evNull->drawAvgSpectrogram(major, minor, 0,ch,evNull->scopeNEvents, nfft, nOverlap, padTo, window, log);
    avgSpecN->SetTitle("Full CH"+TString::Itoa(ch, 10));
    //    avgSpecN->GetZaxis()->SetRangeUser(zmin, zmax);
    TUtil::ranges(avgSpecN, 0, 90, 0, 3);
    cursors[0]->Draw();
    cursors[1]->Draw();
    can->Draw();
    return 1;
}

TNtuple * TAnalyze::integrateAllWithSideband(int channel, int dataset, int major, int minor, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, double sbxmin, double sbxmax, double sbymin, double sbymax, int norm){
  auto ev=new T576Event(50, dataset);
  ev->loadScopeEvent(major, minor, 0);
  int number=0;
  TNtuple *tup=new TNtuple("tup", "tup", "sig:sb");
  //  auto graphs=vector<TGraph*>();
  //if(scopeOrSurf==0){
    //    loadScopeEvent(major, minor, 0);
    for(int i=0;i<ev->scopeNEvents;i++){
      ev->loadScopeEvent(major, minor, i);
      if(ev->isGood==1){
	auto spec = TUtil::FFT::spectrogram(ev->scope->ch[channel], nfft, overlap, zeroPadLength, window, dbFlag);
	if(norm==1){
	  TUtil::normalize(spec, 0, ev->analogBandwidth);
	}
	//spec->Scale(scale);
	spec->SetDirectory(0);
	auto sig=TUtil::integrate(spec, xmin, xmax, ymin, ymax);
	auto sb=TUtil::integrate(spec, sbxmin, sbxmax, sbymin, sbymax);
	tup->Fill(sig, sb);
	delete spec;
	number++;
      }
      //      if(number>=num)break;
    
  }
  return tup;
}
vector<vector<double>> TAnalyze::sidebandSubtractAll(int channel, int dataset,int major, int minor, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, int norm){
  auto ev=new T576Event(50, dataset);
  ev->loadScopeEvent(major, minor, 0);
  int number=0;
  //  TNtuple *tup=new TNtuple("tup", "tup", "sig:err");
  auto graphs=vector<TGraph*>();
  auto tup=vector<vector<double>>(ev->scopeNEvents, vector<double>(2));
    //    loadScopeEvent(major, minor, 0);
    for(int i=0;i<ev->scopeNEvents;i++){
      ev->loadScopeEvent(major, minor, i);
      if(ev->isGood==1){
	auto spec = TUtil::FFT::spectrogram(ev->scope->ch[channel], nfft, overlap, zeroPadLength, window, dbFlag);
	double err=0.;
	if(norm==1){
	  TUtil::normalize(spec, 0, ev->analogBandwidth);
	} //	spec->Scale(scale);
	auto sig=TUtil::sidebandSubtraction2DWithErrors(spec, xmin, xmax, ymin, ymax, err);
	//	auto sb=TUtil::integrate(spec, sbxmin, sbxmax, sbymin, sbymax);
	delete spec;
	//	tup->Fill(sig, err);
	tup[i][0]=sig;
	tup[i][1]=err;
	number++;
      }
      //      if(number>=num)break;
      
  }
  return tup;
}

vector<vector<double>> TAnalyze::sidebandSubtractAllXAxis(int channel, int dataset,int major, int minor, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, int norm){
  auto ev=new T576Event(50, dataset);
  ev->loadScopeEvent(major, minor, 0);
  int number=0;
  //  TNtuple *tup=new TNtuple("tup", "tup", "sig:err");
  auto graphs=vector<TGraph*>();
  auto tup=vector<vector<double>>(ev->scopeNEvents, vector<double>(2));
    //    loadScopeEvent(major, minor, 0);
    for(int i=0;i<ev->scopeNEvents;i++){
      ev->loadScopeEvent(major, minor, i);
      if(ev->isGood==1){
	auto spec = TUtil::FFT::spectrogram(ev->scope->ch[channel], nfft, overlap, zeroPadLength, window, dbFlag);
	double err=0.;
	if(norm==1){
	  TUtil::normalize(spec, 0, ev->analogBandwidth);
	} //	spec->Scale(scale);
	auto sig=TUtil::sidebandSubtractionXAxisWithErrors(spec, xmin, xmax, ymin, ymax, err);
	//	auto sb=TUtil::integrate(spec, sbxmin, sbxmax, sbymin, sbymax);
	delete spec;
	//	tup->Fill(sig, err);
	tup[i][0]=sig;
	tup[i][1]=err;
	number++;
      }
      //      if(number>=num)break;
      
  }
  return tup;
}

TNtuple * TAnalyze::sidebandSubtractAllTuple(int channel, int dataset,int major, int minor, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, int norm){
  auto ev=new T576Event(50, dataset);
  ev->loadScopeEvent(major, minor, 0);
  int number=0;
  TNtuple *tup=new TNtuple("tup", "tup", "sig:err");
  auto graphs=vector<TGraph*>();

    //    loadScopeEvent(major, minor, 0);
    for(int i=0;i<ev->scopeNEvents;i++){
      ev->loadScopeEvent(major, minor, i);
      if(ev->isGood==1){
	auto spec = TUtil::FFT::spectrogram(ev->scope->ch[channel], nfft, overlap, zeroPadLength, window, dbFlag);
	double err=0.;
	if(norm==1){
	  TUtil::normalize(spec, 0, ev->analogBandwidth);
	} //	spec->Scale(scale);
	auto sig=TUtil::sidebandSubtraction2DWithErrors(spec, xmin, xmax, ymin, ymax, err);
	//	auto sb=TUtil::integrate(spec, sbxmin, sbxmax, sbymin, sbymax);
	delete spec;
	tup->Fill(sig, err);
	number++;
      }
      //      if(number>=num)break;
      
  }
  return tup;
}



TH2D * TAnalyze::nullSubtractedSpectrogram(int ch, int major, int minor, int event, int nfft, int nOverlap, int padTo, int window, int log){

  auto ev=new T576Event(50, 1);
  auto evNull=new T576Event(50, 3);

  ev->loadScopeEvent(major, minor, event);
  evNull->loadScopeEvent(major, minor, event);

  auto N =ev->scope->ch[ch]->GetN()<evNull->scope->ch[ch]->GetN()?ev->scope->ch[ch]->GetN():evNull->scope->ch[ch]->GetN();
  auto dGr=TUtil::getNSamplesFrom(ev->scope->ch[ch], 0, N, 0);
  auto dGrN=TUtil::getNSamplesFrom(evNull->scope->ch[ch], 0, N, 0);
 

  
  auto specD=TUtil::FFT::spectrogram(dGr, nfft, nOverlap, padTo, window, log);
  auto specN=TUtil::FFT::spectrogram(dGrN, nfft, nOverlap, padTo, window, log);
  
  auto normDat=TUtil::integrate(specD, specD->GetXaxis()->GetXmin(), specD->GetXaxis()->GetXmin(), specD->GetYaxis()->GetXmin(), specD->GetYaxis()->GetXmin());
  auto normNull=TUtil::integrate(specN, specN->GetXaxis()->GetXmin(), specN->GetXaxis()->GetXmin(), specN->GetYaxis()->GetXmin(), specN->GetYaxis()->GetXmin());

  specD->Scale(1./normDat);
  specN->Scale(1./normNull);
  specD->Add(specN, -1);
  specD->Scale(normDat);
  delete dGr;
  delete dGrN;
  delete ev;
  delete evNull;
  delete specN;
  return specD;
  

  
}


// TH2D * TAnalyze::avgNullSubtractedSpectrogram(int ch, int major, int minor, int nfft, int nOverlap, int padTo, int window, int log){
//   auto ev=new T576Event(50, 1);
//   ev->loadScopeEvent(major, minor, 0);
//   auto avgg=nullSubtractedSpectrogram(ch, major, minor, 0,nfft, nOverlap, padTo, window, log);
//   int num=0;
//   //  for(int i=0;i<ev->scopeNEvents;i++){
//   for(int i=0;i<1;i++){
//     auto temp=nullSubtractedSpectrogram(ch, major, minor, i,nfft, nOverlap, padTo, window, log);
//     avgg->Add(temp);
//     delete temp;
//     num++;
//   }
//   avgg->Scale(1./(double)num);
//   delete ev;
//   return avgg;  
// }


double TAnalyze::getSignalTOA(int ch, int major, int minor, int dataset){
  auto ev=new T576Event(50, dataset);
  ev->loadScopeEvent(major, minor, 0);
  auto t_0=TUtil::getFirstThresholdCrossing(ev->scope->ch[3], .01);
  //distance from beam pipe to center of target. assume that the
  //charges travel at c through the material.
  auto t_1=TUtil::timeOfFlight(ev->beamPipeExit, ev->targetCenter);
  auto compTheta=abs((TUtil::pi/2)-ev->scope->pos[ch].Theta());
  //time of RF propagation through material
  //get the hypotneuse of triangle in direction of tx.
  auto d_2=sqrt((.6*.6)+((.6*tan(compTheta))*(.6*tan(compTheta))));
  //multiply by index of refraction. divide by c.
  //TODO actually solve for the refracted path
  auto t_2=d_2*1.51/TUtil::c_light;

  auto t_3=((ev->targetCenter-ev->scope->pos[ch]).Mag()-d_2)/TUtil::c_light;
  return t_0+t_1+t_2+t_3;


}

double TAnalyze::getSignalDTOA(int ch, int major, int minor, int dataset){
  auto ev=new T576Event(50, dataset);
  ev->loadScopeEvent(major, minor, 0);
  auto t_0=TUtil::timeOfFlight(ev->beamPipeExit, ev->scope->pos[ch]);
  //distance from beam pipe to center of target. assume that the
  //charges travel at c through the material.
  auto t_1=TUtil::timeOfFlight(ev->beamPipeExit, ev->targetCenter);
  auto compTheta=abs((TUtil::pi/2)-ev->scope->pos[ch].Theta());
  //time of RF propagation through material
  //get the hypotneuse of triangle in direction of tx.
  auto d_2=sqrt((.6*.6)+((.6*tan(compTheta))*(.6*tan(compTheta))));
  //multiply by index of refraction. divide by c.
  //TODO actually solve for the refracted path
  auto t_2=d_2*1.51/TUtil::c_light;

  auto t_3=((ev->targetCenter-ev->scope->pos[ch]).Mag()-d_2)/TUtil::c_light;
  return t_1+t_2+t_3-t_0;
}

TH1F* TAnalyze::bootstrapSidebandSubtract(int major, int minor, int channel, int dataset, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, int n, int nbins, int binlow, int binhigh){
  auto vec=TAnalyze::sidebandSubtractAll(channel, dataset, major, minor, nfft, overlap, zeroPadLength, window, dbFlag, xmin, xmax, ymin, ymax, 0);

  auto rann=new TRandom3();
    rann->SetSeed();
    auto histt=new TH1F("asdf", "asdf", nbins, binlow, binhigh);
    for(int i=0;i<n;i++){
      auto avg=0.;
      for(int j=0;j<100;j++){
        auto val=rann->Integer(100);
        avg+=vec[val][0]*1000.;
      }
      histt->Fill(avg/100);
    }
    histt->SetDirectory(0);
    return histt;

  }


TH1F* TAnalyze::bootstrapSidebandSubtractXAxis(int major, int minor, int channel, int dataset, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, int n, int nbins, int binlow, int binhigh){
  auto vec=TAnalyze::sidebandSubtractAllXAxis(channel, dataset, major, minor, nfft, overlap, zeroPadLength, window, dbFlag, xmin, xmax, ymin, ymax, 0);

  auto rann=new TRandom3();
    rann->SetSeed();
    auto histt=new TH1F("asdf", "asdf", nbins, binlow, binhigh);
    for(int i=0;i<n;i++){
      auto avg=0.;
      for(int j=0;j<100;j++){
        auto val=rann->Integer(100);
        avg+=vec[val][0]*1000.;
      }
      histt->Fill(avg/100);
    }
    histt->SetDirectory(0);
    return histt;

  }
