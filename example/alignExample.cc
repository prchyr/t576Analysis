#include "t576/T576Event.hh"
#include "t576/TUtil.hh"

int thing(){
  auto ev=new T576Event();
  //  ev->setInterpGsNs(100);
  vector<TGraph*> graphs;
  for(int i=0;i<20;i++){
    ev->loadScopeEvent(3, 14, i);
    
    graphs.push_back((TGraph*)ev->scope->gr[0]->Clone());

  }

  
  auto alignedGraphs=TUtil::alignMultiple(graphs, 5, 480, 500);

  TUtil::setCoolPalette();
  alignedGraphs[0]->Draw("al PLC");
  for(int i=1;i<alignedGraphs.size();i++){
    alignedGraphs[i]->Draw("l same PLC");
  }
  
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
