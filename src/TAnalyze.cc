/*
  copyright s. prohira and the T576 collaboration 2018
  released under the GNU General Public License version 3
*/
#include <t576/TAnalyze.hh>

int TAnalyze::drawAvgRealNullWithGeom(int ch, int major, int minor, int nfft, int nOverlap, int padTo, int window, int log){
    auto can=new TCanvas("", "", 1200, 400);
    can->Divide(3, 0);
    can->cd(1)->SetRightMargin(.15);

    auto ev=new T576Event(50, 1);
    auto evNull=new T576Event(50, 3);
    //    TUtil::setColdPalette();
    ev->scope->ch[0]->Draw("al PLC");

    auto avgSpec=ev->drawAvgSpectrogram(major, minor, 0,ch,ev->scopeNEvents, nfft, nOverlap, padTo, window, log);
    avgSpec->SetTitle("Data CH"+TString::Itoa(ch, 10));

    TUtil::ranges(avgSpec, 20, 90, 0, 3);

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
    avgSpecN->SetTitle("Null CH"+TString::Itoa(ch, 10));
    avgSpecN->GetZaxis()->SetRangeUser(zmin, zmax);
    TUtil::ranges(avgSpecN, 20, 90, 0, 3);
    cursors[0]->Draw();
    cursors[1]->Draw();
    can->Draw();
    return 1;
}


int TAnalyze::drawAvgRealNull(int ch, int major, int minor, int nfft, int nOverlap, int padTo, int window, int log, int norm){
    auto can=new TCanvas("", "", 1200, 600);
    can->Divide(2, 0);
    can->cd(1)->SetRightMargin(.15);

    auto ev=new T576Event(50, 1);
    auto evNull=new T576Event(50, 3);
    //TUtil::setColdPalette();
    ev->scope->ch[0]->Draw("al PLC");

    auto avgSpec=ev->drawAvgSpectrogram(major, minor, 0,ch,ev->scopeNEvents, nfft, nOverlap, padTo, window, log);
    avgSpec->SetTitle("Data CH"+TString::Itoa(ch, 10));

    TUtil::ranges(avgSpec, 20, 90, 0, 3);

    auto zmax=avgSpec->GetMaximum();
    auto zmin=avgSpec->GetMinimum();
    can->Draw();
    auto cursors=TUtil::drawPeakCursorXY(avgSpec, kRed);
    auto minx=38;
    auto maxx=50.;
    auto miny=1.9;
    auto maxy=2.25;

    can->cd(2)->SetRightMargin(.15);
    auto avgSpecN=evNull->drawAvgSpectrogram(major, minor, 0,ch,evNull->scopeNEvents, nfft, nOverlap, padTo, window, log);

    if(norm==1){
      // auto normDat=TUtil::integrate(avgSpec, avgSpec->GetXaxis()->GetXmin(), avgSpec->GetXaxis()->GetXmin(), avgSpec->GetYaxis()->GetXmin(), avgSpec->GetYaxis()->GetXmin());
      // auto normNull=TUtil::integrate(avgSpecN, avgSpecN->GetXaxis()->GetXmin(), avgSpecN->GetXaxis()->GetXmin(), avgSpecN->GetYaxis()->GetXmin(), avgSpecN->GetYaxis()->GetXmin());
      
      //avgSpecN->Scale(normDat/normNull);
      auto normD=TUtil::norm(avgSpec);
      auto normN=TUtil::norm(avgSpecN);
      avgSpecN->Scale(normD/normN);
    }
    avgSpecN->SetTitle("Null CH"+TString::Itoa(ch, 10));
    avgSpecN->GetZaxis()->SetRangeUser(zmin, zmax);
    TUtil::ranges(avgSpecN, 20, 90, 0, 3);
    cursors[0]->Draw();
    cursors[1]->Draw();
    can->Draw();
    return 1;
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
    avgSpec->SetTitle("Data CH"+TString::Itoa(ch, 10));

    TUtil::ranges(avgSpec, 20, 90, 0, 3);

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
    avgSpecN->SetTitle("Null CH"+TString::Itoa(ch, 10));
    //avgSpecN->GetZaxis()->SetRangeUser(zmin, zmax);
    TUtil::ranges(avgSpecN, 20, 90, 0, 3);
    cursors[0]->Draw();
    cursors[1]->Draw();
    can->Draw();
    return 1;
}


int TAnalyze::drawAvgRealFull(int ch, int major, int minor, int nfft, int nOverlap, int padTo, int window, int log){
    auto can=new TCanvas("", "", 1200, 600);
    can->Divide(2, 0);
    can->cd(1)->SetRightMargin(.15);

    auto ev=new T576Event(50, 1);
    auto evNull=new T576Event(50, 2);
    //TUtil::setColdPalette();
    ev->scope->ch[0]->Draw("al PLC");

    auto avgSpec=ev->drawAvgSpectrogram(major, minor, 0,ch,ev->scopeNEvents, nfft, nOverlap, padTo, window, log);
    avgSpec->SetTitle("Data CH"+TString::Itoa(ch, 10));

    TUtil::ranges(avgSpec, 20, 90, 0, 3);

    auto zmax=avgSpec->GetMaximum();
    auto zmin=avgSpec->GetMinimum();
    can->Draw();
    auto cursors=TUtil::drawPeakCursorXY(avgSpec, kRed);
    auto minx=38;
    auto maxx=50.;
    auto miny=1.9;
    auto maxy=2.25;

    can->cd(2)->SetRightMargin(.15);
    auto avgSpecN=evNull->drawAvgSpectrogram(major, minor, 0,ch,evNull->scopeNEvents, nfft, nOverlap, padTo, window, log);
    avgSpecN->SetTitle("Null CH"+TString::Itoa(ch, 10));
    //    avgSpecN->GetZaxis()->SetRangeUser(zmin, zmax);
    TUtil::ranges(avgSpecN, 20, 90, 0, 3);
    cursors[0]->Draw();
    cursors[1]->Draw();
    can->Draw();
    return 1;
}

TNtuple * TAnalyze::integrateAllWithSideband(int major, int minor, int scopeOrSurf, int channel, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, double sbxmin, double sbxmax, double sbymin, double sbymax, int norm){
  auto ev=new T576Event(50, 0);
  ev->loadScopeEvent(major, minor, 0);
  int number=0;
  TNtuple *tup=new TNtuple("tup", "tup", "sig:sb");
  //  auto graphs=vector<TGraph*>();
  if(scopeOrSurf==0){
    //    loadScopeEvent(major, minor, 0);
    for(int i=0;i<ev->scopeNEvents;i++){
      ev->loadScopeEvent(major, minor, i);
      if(ev->isGood==1){
	auto spec = TUtil::FFT::spectrogram(ev->scope->ch[channel], nfft, overlap, zeroPadLength, window, dbFlag);
	if(norm==1){
	  TUtil::normalize(spec);
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
  }
  return tup;
}

TNtuple * TAnalyze::sidebandSubtractAll(int major, int minor, int scopeOrSurf, int channel, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, int norm){
  auto ev=new T576Event(50, 1);
  ev->loadScopeEvent(major, minor, 0);
  int number=0;
  TNtuple *tup=new TNtuple("tup", "tup", "sig:err");
  auto graphs=vector<TGraph*>();
  if(scopeOrSurf==0){
    //    loadScopeEvent(major, minor, 0);
    for(int i=0;i<ev->scopeNEvents;i++){
      ev->loadScopeEvent(major, minor, i);
      if(ev->isGood==1){
	auto spec = TUtil::FFT::spectrogram(ev->scope->ch[channel], nfft, overlap, zeroPadLength, window, dbFlag);
	double err=0.;
	if(norm=1){
	  TUtil::normalize(spec);
	} //	spec->Scale(scale);
	auto sig=TUtil::sidebandSubtraction2DWithErrors(spec, xmin, xmax, ymin, ymax, err);
	//	auto sb=TUtil::integrate(spec, sbxmin, sbxmax, sbymin, sbymax);
	delete spec;
	tup->Fill(sig, err);
	number++;
      }
      //      if(number>=num)break;
    }
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
