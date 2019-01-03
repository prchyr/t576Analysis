#include "t576/T576Event.hh"
#include "t576/TUtil.hh"

/*
in this example, we take a selection of events from a run and make
a basis out of them. we then take a different event form the same 
run and expand it in the basis, using diffferent numbers of patterns.

we then plot the result, and histogram the difference between the 
reconstucted values and the actual values. 
 */

int thing(){
  auto ev=new T576Event();
  
  vector<TGraph*> graphs;
  for(int i=0;i<20;i++){
    ev->loadScopeEvent(3, 10, i);
    //get a small chunk of the graph to perform SVD on.
    //too long a trace and it'll eat all the memory on your machine (trust me)
    graphs.push_back(TUtil::getChunkOfGraphFine(ev->scope->gr[0], 460, 520));

  }

  //make a matrix from the events, each event will be a row 
  auto M=TUtil::SVD::eventMatrix(graphs);
  //build a basis of those events
  auto basis=TUtil::SVD::buildBasis(M);

  //load a test event to expand in the basis
  ev->loadScopeEvent(3, 10, 15);
  double tmin=460.;
  double tmax=520.;
  TGraph * testGr=TUtil::getChunkOfGraphFine(ev->scope->gr[0], tmin, tmax);

  //expand the event in the basis, using different numbers of patterns.
  TGraph * gOut1= TUtil::SVD::expandInBasis(testGr, basis, 1);
  TGraph * gOut3= TUtil::SVD::expandInBasis(testGr, basis, 3);
  TGraph * gOut10= TUtil::SVD::expandInBasis(testGr, basis, 10);


  gOut1->SetName("1");
  gOut3->SetName("3");
  gOut10->SetName("10");

  //all plotting things now. 

  TCanvas * can=new TCanvas("canvas", "canvas", 1200, 600);
  can->Divide(2, 0);
  can->cd(1);
  TUtil::setCoolPalette();

  testGr->SetTitle("");
  testGr->SetName ("real event");
  testGr->GetXaxis()->SetTitle("time (ns)");
  testGr->GetXaxis()->SetRangeUser(tmin, tmax);
  testGr->SetLineColor(kRed);
  testGr->Draw("al");


  gOut1->Draw("l same PLC");
  gOut3->Draw("l same PLC");
  gOut10->Draw("l same PLC");
  can->cd(1)->BuildLegend(.6, .7, .88, .88, "number of patterns", "l");
  gStyle->SetOptStat(0);
  can->cd(2)->SetGrid();
  TUtil::plotResiduals(gOut1, testGr, 100, -.4, .4)->Draw("l  PLC");
  TUtil::plotResiduals(gOut3, testGr, 100, -.4, .4)->Draw("l same PLC");
  TUtil::plotResiduals(gOut10, testGr, 100, -.4, .4)->Draw("l same PLC");
  auto hist =((TH1F*)can->cd(2)->GetListOfPrimitives()->At(0));
  hist->GetYaxis()->SetRangeUser(0, 475);
  hist->GetXaxis()->SetTitle("#Delta(V_{reconstucted}, V_{actual})");
  hist->SetTitle("Residuals");
  can->Draw();
  return 1;
}


#ifndef __CINT__
void StandaloneApplication(int argc, char** argv) {
  //command to make root run not in root
  thing();
  //  exit(0);
}

//main method-just some ROOT calls to make the app run.
int main(int argc, char** argv) {
  TApplication app("ROOT Application", &argc, argv);
  StandaloneApplication(app.Argc(), app.Argv());
  app.Run();
  return 0;

}
#endif
