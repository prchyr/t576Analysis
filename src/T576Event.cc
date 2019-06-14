/*
copyright s. prohira and the T576 collaboration 2018
released under the GNU General Public License version 3
*/

#include "T576Event.hh"

//using namespace CLHEP;
using namespace std;




//on construction, check that the data directory structure is good and indexed. if not, build an index file for quick grabbing of events.
int T576Event::checkStatus(){
  //where to look for things.
  fInstallDir=getenv("T576_INSTALL_DIR");
  if(fInstallDir==""){
    fInstallDir="/usr/local";
  }
  //find/make the event index file.
  TFile* indexFile;
  indexFile=TFile::Open(fInstallDir+"/share/t576/eventIndex.root");
  if(!indexFile){
    cout<<endl<<"but that's OK! event index is not built yet. it will be built now (it will only happen this one time). you may see zombies below. this is also OK."<<endl<<"building..."<<endl;
    int build=buildEventIndex();
    if(build==0){
      cout<<"can't build event index. set T576_DATA_DIR to your data directory. it is currently set to: "<<getenv("T576_DATA_DIR")<<endl;
      return 0;
    }
    indexFile=TFile::Open(fInstallDir+"/share/t576/eventIndex.root");
  }
  //build the indices for fast lookup.
  fIndexTree=(TTree*)indexFile->Get("indexTree");
  fIndexTree->SetBranchAddress("scopeFilename", &scopeFilename);
  fIndexTree->SetBranchAddress("surfFilename", &surfFilename);
  fIndexTree->SetBranchAddress("major", &major);
  fIndexTree->SetBranchAddress("minor", &minor);
  //fIndexTree->SetBranchAddress("backmajor", &backmajor);
  //fIndexTree->SetBranchAddress("backminor", &backminor);
  fIndexTree->SetBranchAddress("surfNEvents", &surfNEvents);
  fIndexTree->SetBranchAddress("scopeNEvents", &scopeNEvents);
  fIndexTree->SetBranchAddress("subEvNo", &subEvNo);
  fIndexTree->SetBranchAddress("scopeEvNo", &scopeEvNo);
  fIndexTree->SetBranchAddress("surfEvNo", &surfEvNo);
  fIndexTree->SetBranchAddress("scopeTime", &scopeTime);
  fIndexTree->SetBranchAddress("surfTime", &surfTime);
  
  fIndexBuilt=1;
  fIndexTree->BuildIndex("scopeEvNo", "scopeEvNo");
  fScopeEvNoIndex=(TTreeIndex*)fIndexTree->GetTreeIndex();
  fIndexTree->BuildIndex("surfEvNo", "surfEvNo");
  fSurfEvNoIndex=(TTreeIndex*)fIndexTree->GetTreeIndex();
  fIndexTree->BuildIndex("surfTime", "surfTime");
  fSurfTimeIndex=(TTreeIndex*)fIndexTree->GetTreeIndex();
  fIndexTree->BuildIndex("scopeTime", "scopeTime");
  fScopeTimeIndex=(TTreeIndex*)fIndexTree->GetTreeIndex();
  fIndexTree->BuildIndex("major", "minor");
  fMajorMinorIndex=(TTreeIndex*)fIndexTree->GetTreeIndex();

  fIndexTree->GetEntry(fIndexTree->GetEntries()-1);
  fNEntriesScope=scopeEvNo;
  fNEntriesSurf=surfEvNo;

  //load the position/frequency information.

  fRunLog=TFile::Open(fInstallDir+"/share/t576/runLog.root");
  fRunLogTree=(TTree*)fRunLog->Get("runLogTree");
  fRunLogTree->SetBranchAddress("major", &major);
  fRunLogTree->SetBranchAddress("minor", &minor);
  fRunLogTree->SetBranchAddress("backmajor", &backmajor);
  fRunLogTree->SetBranchAddress("backminor", &backminor);
  fRunLogTree->SetBranchAddress("frequency", &frequency);
  fRunLogTree->SetBranchAddress("power", &power);
  fRunLogTree->SetBranchAddress("txOn", &txOn);
  fRunLogTree->SetBranchAddress("polarization", &polarization);
  fRunLogTree->SetBranchAddress("txDist", &txDist);
  fRunLogTree->SetBranchAddress("scopeDist", scope->dist);
  fRunLogTree->SetBranchAddress("txAng", &txAng);
  fRunLogTree->SetBranchAddress("scopeAng", scope->ang);
  fRunLogTree->SetBranchAddress("antennaType", antennaType);
  fRunLogTree->SetBranchAddress("isGood", &isGood);
  fRunLogTree->SetBranchAddress("surfOn", &surfOn);
  fRunLogTree->SetBranchAddress("surfAng", surf->ang);
  fRunLogTree->SetBranchAddress("surfDist", surf->dist);
  //build an index for this tree. we keep the tree open, but save the index too. 
  fRunLogTree->BuildIndex("major", "minor");
  fRunLogIndex=(TTreeIndex*)fRunLogTree->GetTreeIndex();

  fSurfData=(short*) calloc(60000000, sizeof(short));
  fSurfTimes=(double*) calloc(20000, sizeof(double));

  beamPipeExit.SetXYZ(.6, 0., -3.6);
  targetFront.SetXYZ(.6, 0., -2.);
  targetBack.SetXYZ(.6, 0., 2.);
  beamDump.SetXYZ(.6, 0., 4.);

}











//load a scope event using a run major/minor combination and event within that file

int T576Event::loadScopeEvent(int run_major, int run_minor,int event, bool remove_dc_offset){
  if(fIndexBuilt==0){
    cout<<"event index not built yet. building..."<<endl;
    buildEventIndex();
  }
  TString top_dir = getenv("T576_DATA_DIR");

  if(fUSE_FILTERED_DATA==1){
    top_dir=getenv("T576_FILTERED_DATA_DIR");
    if(top_dir==""){
      cout<<"T576_FILTERED_DATA_DIR not set. please set this flag so that the filtered data can be found. "<<endl;
      return (0);
    }
  }
  
  if(fUSE_FILTERED_DATA==2){
    top_dir=getenv("T576_NULL_DATA_DIR");
    if(top_dir==""){
      cout<<"T576_FILTERED_DATA_DIR not set. please set this flag so that the filtered data can be found. "<<endl;
      return (0);
    }
  }
 
 if(fUSE_FILTERED_DATA==3){
    top_dir=getenv("T576_FILTERED_NULL_DATA_DIR");
    if(top_dir==""){
      cout<<"T576_FILTERED_DATA_DIR not set. please set this flag so that the filtered data can be found. "<<endl;
      return (0);
    }
  }
  

  
  if(top_dir==""){
   cout<<"T576_DATA_DIR not set. please set this flag so that the data can be found. this should be the top directory inside of which is py/ and root/."<<endl;
    return (0);
  }
  
  TString directory=top_dir+"/root/";

  
  fIndexTree->SetTreeIndex(fMajorMinorIndex);
  fIndexTree->GetEntry(fIndexTree->GetEntryNumberWithBestIndex(run_major, run_minor));
  fRunLogTree->GetEntry(fRunLogTree->GetEntryNumberWithBestIndex(run_major, run_minor));

  txPos.SetMagThetaPhi(txDist, txAng*pi/180., pi);

  for(int i=0;i<4;i++){
    scope->pos[i].SetMagThetaPhi(scope->dist[i], scope->ang[i]*pi/180., pi);
  }

  TString thisScopeFilename=scopeFilename->Data();
  //if the file is already open, save some time
  if(thisScopeFilename!=fScopeFilename){

    fScopeFilename=thisScopeFilename;

    //open the file
    fEventFile->Close();
    fEventFile=TFile::Open(directory+thisScopeFilename);
    fEventTree=(TTree*)fEventFile->Get("tree");
  
    //set addresses
  
    fEventTree->SetBranchAddress("ch1", scope->dat[0]);
    fEventTree->SetBranchAddress("ch2", scope->dat[1]);
    fEventTree->SetBranchAddress("ch3", scope->dat[2]);
    fEventTree->SetBranchAddress("ch4", scope->dat[3]);
    fEventTree->SetBranchAddress("time", scope->time);
    fEventTree->SetBranchAddress("timestamp", &timestamp);

  
  }
  
  if(major!=run_major&&minor!=run_minor){
    cout<<"major/minor combination doesn't exist"<<endl;
    return 0;
  }
  if(loadScopeEvent(scopeEvNo+event, remove_dc_offset)==1){
    getCharge(scope->ch[3]);
    //fEventFile->Close();
    return 1;
  }
  else return 0;
  
  // if(event>fEventTree->GetEntries()){
  //   cout<<"event number too high! this major/minor combination only contains "<<fEventTree->GetEntries()<<" events."<<endl<<"select another event number, or call loadScopeEvent(event) or loadSurfEvent(event) without major/minor to use overall event index"<<endl;
  //   return (-1);
  // }

  // fEventTree->GetEntry(event);

  //  //check the length of the record.
  // auto length=sizeof(scope->time)/sizeof(*scope->time);

  // //fill the event graphs for the scope->
  // //fix the first and last values, which were recorded incorrectly
  // scope->time[0]=0.;
  // scope->time[19999]=scope->time[19998]+.05;
  // for(int i=0;i<4;i++){
  //   TGraph * graph=new TGraph(length, scope->time, scope->dat[i]);

  //   if(fInterpGSs>0.){
  //     getInterpolatedGraph(graph, scope->ch[i]);
  //   }
  //   else{
  //     *scope->ch[i]=*graph;
  //   }
  //   delete(graph);
  // }
  
  // //delete(tree);
  // fEventFile->Close();
  // //delete(file);
  // getCharge(scope->ch[3]);
  
  // //  delete (files);
  
  
  // return 1;
  }


//load an event from the global scope event number index 
int T576Event::loadScopeEvent(int event, bool remove_dc_offset){

  if(fIndexBuilt==0){
    cout<<"event index not built yet. building..."<<endl;
    buildEventIndex();
  }
  
  TString top_dir = getenv("T576_DATA_DIR");

  if(fUSE_FILTERED_DATA==1){
    top_dir=getenv("T576_FILTERED_DATA_DIR");
    if(top_dir==""){
      cout<<"T576_FILTERED_DATA_DIR not set. please set this flag so that the filtered data can be found. "<<endl;
      return (0);
    }
  }
  
  if(fUSE_FILTERED_DATA==2){
    top_dir=getenv("T576_NULL_DATA_DIR");
    if(top_dir==""){
      cout<<"T576_FILTERED_DATA_DIR not set. please set this flag so that the filtered data can be found. "<<endl;
      return (0);
    }
  }
 
 if(fUSE_FILTERED_DATA==3){
    top_dir=getenv("T576_FILTERED_NULL_DATA_DIR");
    if(top_dir==""){
      cout<<"T576_FILTERED_DATA_DIR not set. please set this flag so that the filtered data can be found. "<<endl;
      return (0);
    }
  }
  

  
  if(top_dir==""){
   cout<<"T576_DATA_DIR not set. please set this flag so that the data can be found. this should be the top directory inside of which is py/ and root/."<<endl;
    return (0);
  }

  TString directory=top_dir+"/root/";

  
  //cout<<"asdf"<<endl;
  //  fIndexTree->SetBranchAddress("scopeFilename", &scopeFilename);
  fIndexTree->SetTreeIndex(fScopeEvNoIndex);
  fIndexTree->GetEntry(fIndexTree->GetEntryNumberWithBestIndex(event, event));

  //loadScopeEvent(major, minor, subEvNo);
  //return 1;
  //cout<<fIndexBuilt<<endl;
  fRunLogTree->GetEntry(fRunLogTree->GetEntryNumberWithBestIndex(major, minor));
  //cout<<txDist<<" "<<txAng<<endl;
  txPos.SetMagThetaPhi(txDist, txAng*pi/180., pi);
  
  for(int i=0;i<4;i++){
    scope->pos[i].SetMagThetaPhi(scope->dist[i], scope->ang[i]*pi/180., pi);
  }
  if(major<1)return 0;
  //  fIndexTree->GetEntry(0);
  //cout<<major<<" "<<minor<<endl;
  TString thisScopeFilename=scopeFilename->Data();
  //open the file
  //cout<<directory+thisScopeFilename<<endl<<" "<<fScopeFilename<<endl;
  if(thisScopeFilename!=fScopeFilename){

    fEventFile->Close();
    // cout<<thisScopeFilename<<endl;
    fEventFile=TFile::Open(directory+thisScopeFilename);
    fEventTree=(TTree*)fEventFile->Get("tree");

    //  memset(scope->ch, 0, sizeof(scope->ch));
    //memset(scope->gr, 0, sizeof(scope->gr));

    fEventTree->SetBranchAddress("ch1", scope->dat[0]);
    fEventTree->SetBranchAddress("ch2", scope->dat[1]);
    fEventTree->SetBranchAddress("ch3", scope->dat[2]);
    fEventTree->SetBranchAddress("ch4", scope->dat[3]);
    fEventTree->SetBranchAddress("time", scope->time);
    fEventTree->SetBranchAddress("timestamp", &timestamp);

    fScopeFilename=thisScopeFilename;

  }

  if(event>=fNEntriesScope){
    cout<<"event number too high!"<<endl;
    return 0;
  }

  //cout<<"hellpo"<<endl;
  int subEvNo=event-scopeEvNo;
  fEventTree->GetEntry(subEvNo);
  //cout<<"hi"<<endl;
  //check the length of the record.
  auto length=4999;
  if(fUSE_FILTERED_DATA==0){
    length=sizeof(scope->time)/sizeof(*scope->time);
    if(length!=20000)cout<<length;
  }
  //  cout<<length;
  //cout<<"l337"<<endl;
  //cout<<major<<" "<<minor<<" "<<subEvNo<<" "<<length<<" "<<fEventTree->GetEntries()<<" "<<scope->dat[1][18]<<endl;
  //fill the event graphs for the scope->
  //fix the first and last values, which were recorded incorrectly
  if(fUSE_FILTERED_DATA==0){
    scope->time[0]=0.;
    scope->time[19999]=scope->time[19998]+.05;
  }
  for(int i=0;i<4;i++){
    if(fUSE_FILTERED_DATA==0){
      scope->delays[i]=(scope->cableLengths[i]*ft)/(c_light*scope->velocityFactor);
    }
    else{
      scope->delays[i]=0.;
    }
    //cout<<"l352"<<endl;
    TGraph * tempGr=new TGraph(length, scope->time, scope->dat[i]);
    TGraph * graph=TUtil::delayGraph(tempGr, scope->delays[i]);
    //cout<<"l355"<<endl;
    if(fUSE_FILTERED_DATA==0){
      if(remove_dc_offset==true){
	TUtil::removeMeanInPlace(graph, 0., 300.);
      }
      //cout<<"l359"<<endl;
      
      if(fInterpGSs>0.){
	getInterpolatedGraph(graph, scope->ch[i], fInterpGSs);
      }
      else{
	*scope->ch[i]=*graph;
      }
    }
    else{
      //*scope->ch[i]=*tempGr;//cout<<"l370"<<endl;
      getInterpolatedGraph(graph, scope->ch[i], fInterpGSs);
    }
    // cout<<"l372"<<endl;
    scope->ch[i]->SetTitle("");
    scope->ch[i]->SetName("ch"+TString::Itoa(i, 10));
    scope->ch[i]->GetXaxis()->SetTitle("Time (ns)");
    scope->ch[i]->GetYaxis()->SetTitle("Volts (V)");
    scope->ch[i]->GetYaxis()->SetTitleOffset(1.15);
    scope->ch[i]->GetHistogram()->SetName("");
    
   delete(graph);
   delete(tempGr);

  }
  //cout<<"l384"<<endl;
  //  cout<<"here"<<endl;
  //  if(subEvNo==fEventTree->GetEntries())fEventFile->Close();
  //delete(file);
  getCharge(scope->ch[3]);
  
  //  delete (files);

  if(charge<.1){
    isGood=0;
  }

  fScopeLoaded=true;  
  return 1;
}



int T576Event::loadSurfEvent(int run_major, int run_minor, int event, bool remove_dc_offset){

  if(fIndexBuilt==0){
    cout<<"event index not built yet. building..."<<endl;
    buildEventIndex();
  }
  
  TString top_dir = getenv("T576_DATA_DIR");

  if(top_dir==""){
    cout<<"T576_DATA_DIR not set. please set this flag so that the data can be found. this should be the top directory inside of which is py/ and root/."<<endl;
    return (0);
  }

  TString directory=top_dir+"/py/dat/";

  fIndexTree->SetTreeIndex(fMajorMinorIndex);
  fIndexTree->GetEntry(fIndexTree->GetEntryNumberWithBestIndex(run_major, run_minor));


  fRunLogTree->GetEntry(fRunLogTree->GetEntryNumberWithBestIndex(major, minor));

  txPos.SetMagThetaPhi(txDist, txAng*pi/180., pi);
  
  for(int i=0;i<12;i++){
    surf->pos[i].SetMagThetaPhi(surf->dist[i], surf->ang[i]*pi/180., pi);
  }
  if(run_major!=major &&run_minor!=minor){
    cout<<"no surf file for requested major/minor combination."<<endl;
    return 0;
  }
  if(strncmp(surfFilename->Data(), "null", 4)==0){
    cout<<"no surf events for this run."<<endl;
    return 0;
  }
  if(major<1)return 0;
  if(event>surfNEvents){
    cout<<"event number too large for this major/minor combination"<<endl;
    return 0;
  }

  if(loadSurfEvent(surfEvNo+event, remove_dc_offset)==1)return 1;
  else return 0;
//   TString thisSurfFilename=surfFilename->Data();
//   //open the file

//   if(thisSurfFilename!=fSurfFilename){
//     //load the event file
//     TString openf=directory+thisSurfFilename;

//     fDataset = cnpy::npz_load(openf.Data());

//     //load the data. it is stored as an array of shorts of length (N events) x (12 channels ) x (1024 samples)

//     fSurfDataArray = fDataset["data"];

//     fSurfData=fSurfDataArray.data<short>();

//     cnpy::NpyArray times_arr = fDataset["times"];
//     fSurfTimes = times_arr.data<double>();

//     fSurfFilename=thisSurfFilename;
// }


//   if(event>=fNEntriesSurf){
//     cout<<"event number too high!"<<endl;
//     return 0;
//   }


//   subEvNo=event;//-surfEvNo;

//   surfTime=fSurfTimes[subEvNo];

//   //increment and index variables for event number and channel number
//   int ev=12*1024;
//   int len=1024;
//   int index1=subEvNo*ev;

//   //cast the shorts to doubles for plotting and mathing  
//   double data[1024];

//   for(int i=0;i<12;i++){
//     int index2=i*len;//i is channel number
//     copy(fSurfData+index1+index2, fSurfData+index1+index2+len, surf->dat[i]);

//     TGraph * graph=new TGraph(len, surf->time, surf->dat[i]);

//     if(fInterpGSs>0.){
//       getInterpolatedGraph(graph, surf->ch[i]);
//     }
//     else{
//       *surf->ch[i]=*graph;
//     }
//    delete(graph);
//   }


}




//load an event from the global surf event number index 
int T576Event::loadSurfEvent(int event, bool remove_dc_offset){

  if(fIndexBuilt==0){
    cout<<"event index not built yet. building..."<<endl;
    buildEventIndex();
  }
  
  TString top_dir = getenv("T576_DATA_DIR");

  if(top_dir==""){
    cout<<"T576_DATA_DIR not set. please set this flag so that the data can be found. this should be the top directory inside of which is py/ and root/."<<endl;
    return (0);
  }

  TString directory=top_dir+"/py/dat/";

 
  fIndexTree->SetTreeIndex(fSurfEvNoIndex);
  fIndexTree->GetEntry(fIndexTree->GetEntryNumberWithBestIndex(event, event));

  fRunLogTree->GetEntry(fRunLogTree->GetEntryNumberWithBestIndex(major, minor));

  //  while(surfFilename=="null"){
  //fIndexTree->GetEntry(
  txPos.SetMagThetaPhi(txDist, txAng*pi/180., pi);
  
  for(int i=0;i<12;i++){
    surf->pos[i].SetMagThetaPhi(surf->dist[i], surf->ang[i]*deg, pi);
  }
  if(major<1)return 0;

  TString thisSurfFilename=surfFilename->Data();
  //open the file
  //cout<<thisSurfFilename<<endl<<" "<<event<<endl;
  if(thisSurfFilename!=fSurfFilename){
    //load the event file
    TString openf=directory+thisSurfFilename;
    //    fDataset.clear();
    // free(fSurfData);
    // free(fSurfTimes);
    // fSurfData=(short*) calloc(60000000, sizeof(short));
    // fSurfTimes=(double*) calloc(20000, sizeof(double));

    fDataset = cnpy::npz_load(openf.Data());

    //load the data. it is stored as an array of shorts of length (N events) x (12 channels ) x (1024 samples)


    fSurfDataArray = fDataset["data"];
    //auto temp=(short*) calloc(fSurfDataArray.num_bytes(), sizeof(short));
    //auto temp=fSurfDataArray.data<short>();
    fSurfData=fSurfDataArray.data<short>();
    //memcpy(fSurfData, temp, fSurfDataArray.num_bytes());
    //    free(temp);
    //delete temp;
    fTimesArray = fDataset["times"];
    fSurfTimes = fTimesArray.data<double>();

    // delete fDataset;
    // delete fSurfDataArray;
    // delete fTimesArray;
    fSurfFilename=thisSurfFilename;
}


  if(event>=fNEntriesSurf){
    cout<<"event number too high!"<<endl;
    return 0;
  }

  subEvNo=event-surfEvNo;

  surfTime=fSurfTimes[subEvNo];

  //increment and index variables for event number and channel number
  int ev=12*1024;
  int len=1024;
  int index1=subEvNo*ev;

  //cast the shorts to doubles for plotting and mathing  
  double data[1024];

  for(int i=0;i<12;i++){
    int index2=i*len;//i is channel number
    copy(fSurfData+index1+index2, fSurfData+index1+index2+len, surf->dat[i]);
    //feet to meter conversion
    surf->delays[i]=(surf->cableLengths[i]*ft)/(c_light*surf->velocityFactor);
    //mV to V conversion
    for(int j=0;j<len;j++){
      surf->dat[i][j]*=.001;
    }
    TGraph * graph=new TGraph(len, TUtil::makeIndices(len, 1./3.2, surf->delays[i]), surf->dat[i]);

    TGraph *grChunk=(TGraph*)graph->Clone();//TUtil::getChunkOfGraph(graph, 0, 290);
    grChunk->SetName("chunk");
    grChunk->SetTitle("chunk");

    graph->SetName("grrr");
    graph->SetTitle("grrr");
    if(remove_dc_offset==true){
      TUtil::removeMeanInPlace(grChunk, 0., 30.);
    }

    //surf channel mappings are 0, 11,10, 9,...1
    //use sinc interpolation for sampling near nyquist.
    if(fInterpGSs>0.){
      if(i==0){
    	getInterpolatedGraph(grChunk, surf->ch[i], fInterpGSs, 2);
      }
      else{
    	getInterpolatedGraph(grChunk, surf->ch[12-i], fInterpGSs, 2);
      }
    }
    else{
      if(i==0){
    	*surf->ch[i]=*grChunk;
      }
      else{
    	*surf->ch[12-i]=*grChunk;
      }
    }

    // if(fInterpGSs>0.){
    //   getInterpolatedGraph(grChunk, surf->ch[i], fInterpGSs);
    // }
    // else{
    //   *surf->ch[i]=*grChunk;
    // }

    surf->ch[i]->SetTitle("");
    surf->ch[i]->SetName("ch"+TString::Itoa(i, 10));
    surf->ch[i]->GetXaxis()->SetTitle("Time (ns)");
    surf->ch[i]->GetYaxis()->SetTitle("Volts (V)");
    surf->ch[i]->GetYaxis()->SetTitleOffset(1.15);
    surf->ch[i]->GetHistogram()->SetName("");
    
    delete(grChunk);
    delete(graph);
  }

  fSurfLoaded=true;
  //getCharge();  
  return 1;
}


double T576Event::getCharge(TGraph *ict){
  //TUtil::removeMeanInPlace(ict, 0., 400.)
  if(fUSE_FILTERED_DATA==0){
    charge = .4*TUtil::integrate(ict, 435, 550);
  }
  else{
    charge = .4*TUtil::integrate(ict, 5, 90);
  }
  return charge;
}


vector<TGraph *> T576Event::draw(int major, int minor, int scopeOrSurf, int channel, int num, double align, double tLow, double tHigh, TString drawOption){

  int number=0;
  auto graphs=vector<TGraph*>();
  if(scopeOrSurf==0){
    loadScopeEvent(major, minor, 0);
    tHigh=tHigh==999999?tHigh=scope->ch[channel]->GetY()[scope->ch[channel]->GetN()-1]:tHigh;
    for(int i=0;i<scopeNEvents;i++){
      loadScopeEvent(major, minor, i);
      if(isGood){
	graphs.push_back(TUtil::getChunkOfGraphFast(scope->ch[channel], tLow, tHigh));
	number++;
      }
      if(number>=num)break;
    }
    if(align>0.){
      auto aligned=TUtil::alignMultiple(graphs, align);
    //    graphs.clear();
      TUtil::draw(aligned);
      return aligned;
    }
    else{
      return graphs;
    }

  }
  else{
    loadSurfEvent(major, minor, 0);
    tHigh=tHigh==999999?tHigh=surf->ch[channel]->GetY()[surf->ch[channel]->GetN()-1]:tHigh;
    for(int i=0;i<surfNEvents;i++){
      loadSurfEvent(major, minor, i*3);
      if(isGood){
	graphs.push_back(TUtil::getChunkOfGraphFast(surf->ch[channel], tLow, tHigh));
	number++;
      }
      if(number>=num)break;
    }
    if(align>0.){
      auto aligned=TUtil::alignMultiple(graphs, align);
    //    graphs.clear();
      TUtil::draw(aligned);
      return aligned;
    }
    else{
      return graphs;
    }
    
  }


  
}



TGraph * T576Event::drawAvg(int major, int minor, int scopeOrSurf, int channel, int num, double align, double tLow, double tHigh, TString drawOption){
  int number=0;
  auto graphs=vector<TGraph*>();
  if(scopeOrSurf==0){
    for(int i=0;i<scopeNEvents;i++){
      loadScopeEvent(major, minor, i);
      if(isGood){
	graphs.push_back(TUtil::getChunkOfGraphFast(scope->ch[channel], tLow, tHigh));
	number++;
      }
      if(number>=num)break;
    }
    auto avg=TUtil::alignMultipleAndAverage(graphs, align, tLow, tHigh);
    graphs.clear();
    avg->Draw(drawOption);
    return avg;
  }
  else{
    for(int i=0;i<surfNEvents;i++){
      loadSurfEvent(major, minor, i*3);
      if(isGood){
	graphs.push_back((TGraph*)surf->ch[channel]->Clone());
	number++;
      }
      if(number>=num)break;
    }
    auto temp=TUtil::alignMultipleAndTruncate(graphs, align, tLow, tHigh);
    auto avg=TUtil::avgGraph(temp);
    graphs.clear();
    temp.clear();
    avg->Draw(drawOption);
    return avg;
  }


  
}


TH2D * T576Event::drawAvgSpectrogram(int major, int minor, int scopeOrSurf, int channel, int num, Int_t binsize , Int_t overlap, Int_t zero_pad_length, int win_type, int dbFlag){
  int number=0;
  auto specs=vector<TH2D*>();
  if(scopeOrSurf==0){
    for(int i=0;i<scopeNEvents;i++){
      loadScopeEvent(major, minor, i);
      if(isGood){
	//graphs.push_back(TUtil::getChunkOfGraphFast(scope->ch[channel], tLow, tHigh));
	specs.push_back(TUtil::FFT::spectrogram(scope->ch[channel], binsize, overlap, zero_pad_length, win_type, dbFlag));
	number++;
      }
      if(number>=num)break;
    }
    auto avg=TUtil::FFT::avgSpectrograms(specs);
      specs.clear();
      
      avg->Draw("colz");
      return avg;
  }
  // else{
  //   for(int i=0;i<surfNEvents;i++){
  //     loadSurfEvent(major, minor, i*3);
  //     if(isGood){
  // 	graphs.push_back((TGraph*)surf->ch[channel]->Clone());
  // 	number++;
  //     }
  //     if(number>=num)break;
  //   }
  //   auto temp=TUtil::alignMultipleAndTruncate(graphs, align, tLow, tHigh);
  //   auto avg=TUtil::avgGraph(temp);
  //   graphs.clear();
  //   temp.clear();
  //   avg->Draw(drawOption);
  //   return avg;
  // }


  
}


// int T576Event::getInterpolatedGraph(TGraph * inGraph, TGraph *outGraph){
//   ROOT::Math::Interpolator interp(inGraph->GetN(), ROOT::Math::Interpolation::kAKIMA);
//   interp.SetData(inGraph->GetN(), inGraph->GetX(), inGraph->GetY());

//   //get dt, assuming even sampling.
//   double inDt=inGraph->GetX()[50]-inGraph->GetX()[49];
//   double inGsPerS=1./inDt;

//   double outDt=1./fInterpGSs;

//   int samps=(int) (inGraph->GetN()*(fInterpGSs/inGsPerS));

//   vector<double> xx, yy;
//   for(int i=0;i<samps;i++){
//     double time = i*outDt;
//     if(time>inGraph->GetX()[inGraph->GetN()-1])continue;
//     xx.push_back(time);
//     yy.push_back(interp.Eval(time));
//   }

//   auto tempGr=new TGraph(xx.size(), &xx[0], &yy[0]);
//   *outGraph=*tempGr;
//   delete(tempGr);
			 

//   return 1;
// }




int T576Event::buildEventIndex(int force){
  TString install_dir=getenv("T576_INSTALL_DIR");
  if(install_dir==""){
    install_dir="/usr/local";
  }
  

  ifstream indexFile;
  indexFile.open(install_dir+"/share/t576/eventIndex.root");
  if(indexFile&&force==0){
    cout<<"event index already built."<<endl;
    return 1;
  }
  else if(indexFile&&force==1){
    cout<<"event index already built, forcing rebuild.."<<endl;
  }

  //cout<<"asdf"<<endl;
  TFile * index= new TFile(install_dir+"/share/t576/eventIndex.root", "recreate");
  TTree *indexTree = new TTree("indexTree", "index of events");
  TString scopeFilename, surfFilename;
  int maj=0, min=0, subEvNo=0, scopeEvNo=0, surfEvNo=0, scopeNEvents=0, surfNEvents=0;
  ULong64_t tstampScope, tstampSurf;
  indexTree->Branch("scopeFilename", &scopeFilename);
  indexTree->Branch("surfFilename", &surfFilename);
  indexTree->Branch("major", &maj);
  indexTree->Branch("minor", &min);
  indexTree->Branch("scopeNEvents", &scopeNEvents);
  indexTree->Branch("surfNEvents", &surfNEvents);
  indexTree->Branch("subEvNo", &subEvNo);
  indexTree->Branch("scopeEvNo", &scopeEvNo);
  indexTree->Branch("surfEvNo", &surfEvNo);
  indexTree->Branch("scopeTime", &tstampScope);
  indexTree->Branch("surfTime", &tstampSurf);
  indexTree->Branch("tstamp", &tstampScope);

  TString top_dir = getenv("T576_DATA_DIR");
  TString directory=top_dir+"/root/";
  
  TSystemDirectory dir(directory, directory);
  TList *files = dir.GetListOfFiles();
  //    cout<<files->GetEntries()<<endl;
  //find the scopeFilename
  
  if(files){
    files->Sort();
    int nfiles=files->GetEntries();
    TString fname;
    TString ext = ".root";
    //TObject *obj;
    //    TIter next(files);
    for(int i=0;i<nfiles;i++){
      TSystemFile *file=(TSystemFile*)files->At(i);
      //delete(filen);
      fname = file->GetName();

      if(!file->IsDirectory()&&!file->IsZombie()){//is the file there
	TString majorstr=fname(17);//get the run major
	maj=majorstr.Atoi(); 
	//cout<<majorstr<<" ";
	TString substr1=fname(19, 999);
	TString minorstr=substr1(0, substr1.First("."));//get run minor
	//cout<<minorstr<<endl;
	if(minorstr.IsDec()){//check that it isn't 'test' or something
	  min=minorstr.Atoi();
	  //cout<<min<<endl;
	  auto inFile=TFile::Open(directory+fname);//open the file
	  if(inFile){
	    TTree *tree = (TTree*)inFile->Get("tree");
	    if(tree){
	      tree->SetBranchAddress("timestamp", &tstampScope);
	      scopeNEvents=tree->GetEntries();
	      cout.flush()<<fname<<"                  \r";	      
	      //	      for(int j=0;j<nentries;j++){
		//cout<<scopeFilename<<" "<<major<<" "<<minor<<" "<<subEvNo<<" "<<scopeEvNo<<endl;
		tree->GetEntry(0);
		//		subEvNo=j;
		//	      scopeFilename=const_cast<char*>(fname.Data());
		scopeFilename=fname;
		//indexTree->Fill();

		surfFilename=getSurfFilename(maj, min);
		if(surfFilename!="null"){
		  TString sfname=top_dir+"/py/dat/"+surfFilename;
		  cnpy::npz_t dataset = cnpy::npz_load(sfname.Data());
		  //load the data. it is stored as an array of shorts of length (N events) x (12 channels ) x (1024 samples)
		  cnpy::NpyArray data_arr = dataset["data"];
		  short* data_short = data_arr.data<short>();
		  //load the times, stored as double
		  cnpy::NpyArray times_arr = dataset["times"];
		  double * times = times_arr.data<double>();
		  
		  surfNEvents=data_arr.shape[0]/(12*1024);
		  tstampSurf=times[0];

		}
		else{
		  surfNEvents=0;
		  tstampSurf=0;

		}
		indexTree->Fill();
		surfEvNo+=surfNEvents;
		scopeEvNo+=scopeNEvents;
	    }
	  }
	  //inFile->Close();
	
	  delete(inFile);
	  //delete(tree);
	}
      
      }
      delete(file);
    }
  
  }
  //delete (file);
  //}
  else{
    cout<<"no files! check directory."<<endl;
    return 0;
  }

  index->Write();
  index->Close();
  //delete(indexTree);
  delete(index);
  cout.flush();
  cout<<endl<<"index built."<<endl;
  
  fIndexBuilt=1;
  
  return 1;
}

TString T576Event::getSurfFilename(int inMaj, int inMin){
  TString top_dir = getenv("T576_DATA_DIR");
  TString directory=top_dir+"/py/dat/";
  
  TSystemDirectory dir(directory, directory);
  TList *files = dir.GetListOfFiles();
  if(files){
    files->Sort();
    int nfiles=files->GetEntries();
    TString fname;
    TString ext = ".npz";
    //TObject *obj;
    //    TIter next(files);
    for(int i=0;i<nfiles;i++){
      TSystemFile *file=(TSystemFile*)files->At(i);
      //delete(filen);
      fname = file->GetName();

      if(!file->IsDirectory()){//is the file there
	TString majorstr=fname(17);//get the run major
	int maj=majorstr.Atoi();
	if(inMaj!=maj)continue;
	//cout<<majorstr<<" ";
	TString substr1=fname(19, 999);
	TString minorstr=substr1(0, substr1.First("."));//get run minor
	//cout<<minorstr<<endl;
	if(minorstr.IsDec()){//check that it isn't 'test' or something
	  int min=minorstr.Atoi();
	  if(inMin!=min)continue;
	  return fname;
	}
      }
    }
  }
  return "null";
}




/***************************************************************

scope utilities

**************************************************************/

int T576Event::drawGeom(int scopeChan, int surfChan){
  

  auto dummy=new TGraph();

  dummy->SetPoint(dummy->GetN(), -8, -8);
  dummy->SetPoint(dummy->GetN(), 8, -8);
  dummy->SetPoint(dummy->GetN(), 8, 8);
  dummy->SetPoint(dummy->GetN(), -8, 8);
  dummy->SetMarkerSize(0);
  dummy->Draw("ap");
  dummy->GetXaxis()->SetRangeUser(-8, 8);
  dummy->GetYaxis()->SetRangeUser(-8, 8);
  dummy->GetXaxis()->SetTitle("z (m)");
  dummy->GetYaxis()->SetTitle("x (m)");

  double scopex[4], scopey[4], surfx[12], surfy[12];

  for(int i=0;i<12;i++){
    surfx[i]=surf->pos[i].Z();
    surfy[i]=surf->pos[i].X();
    //        cout<<surf[i].x<<" "<<surf[i].y<<" "<<surf[i].z<<endl;
  }


  
  for(int i=0;i<3;i++){
    scopex[i]=scope->pos[i].Z();
    scopey[i]=scope->pos[i].X();
  }

  

  auto txgraph=new TGraph();
  txgraph->SetPoint(0, txPos.Z(), txPos.X());
  txgraph->SetMarkerColor(kBlue);
  txgraph->SetMarkerStyle(21);
  txgraph->SetMarkerSize(1.5);
  txgraph->Draw("p same");

 

  double bdx[5], bdy[5];
  bdx[0] = 2;
  bdx[1]=5.6;
  bdx[2]=5.6;
  bdx[3]=2;
  bdx[4]=2;
  
  bdy[0]=2.44;
  bdy[1]=2.44;
  bdy[2]=-1.24;
  bdy[3]=-1.24;
  bdy[4]=2.44;
  auto bd=new TPolyLine(5, bdx, bdy);
  bd->SetLineColor(kBlack);
  bd->Draw("l same");

  double wallx[5], wally[5];
  wallx[0] = 5.6;
  wallx[1]=5.6;
  
  wally[0]=8;
  wally[1]=-8;

  auto wall=new TPolyLine(2, wallx, wally);
  wall->SetLineColor(kBlack);
  wall->SetLineStyle(3);
  wall->Draw("l same");

  double bd2x[5], bd2y[5];
  bd2x[0] = 3.83;
  bd2x[1]=3.83;
  
  bd2y[0]=-1.24;
  bd2y[1]=2.44;

  auto bd2=new TPolyLine(2, bd2x, bd2y);
  bd2->SetLineColor(kBlack);
  bd2->Draw("l same");

  double blockx[5], blocky[5];
  blockx[0] = -6;
  blockx[1]=-5;
  blockx[2]=-5.5;
  blockx[3]=-6.5;
  blockx[4]=-6;
  
  blocky[0]=-3;
  blocky[1]=-5;
  blocky[2]=-5.5;
  blocky[3]=-3.5;
  blocky[4]=-3;
  auto block=new TPolyLine(5, blockx, blocky);
  block->SetLineColor(kBlack);
  block->Draw("l same");


  double supportx[5], supporty[5];
  supportx[0] = -7.2075;
  supportx[1]=2;
  supportx[2]=2;
  supportx[3]=-7.2075;
  supportx[4]=-7.2075;
  
  supporty[0]=1.517;
  supporty[1]=1.517;
  supporty[2]=-.317;
  supporty[3]=-.317;
  supporty[4]=1.517;
  auto support=new TPolyLine(5, supportx, supporty);
  support->SetLineColor(kBlack);
  support->Draw("l same");

    double bpx[5], bpy[5];
  bpx[0] = -8.;
  bpx[1]=-3.6;
  bpx[2]=-3.6;
  bpx[3]=-8;

  bpy[0]=.5;
  bpy[1]=.5;
  bpy[2]=.7;
  bpy[3]=.7;
  bpy[4]=0;
  auto bp=new TPolyLine(4, bpx, bpy);
  bp->SetLineColor(kBlack);
  bp->SetLineStyle(2);
  bp->Draw("l same");


  double tgtx[5], tgty[5];
  tgtx[0] = -2.;
  tgtx[1]=2.;
  tgtx[2]=2.;
  tgtx[3]=-2.;
  tgtx[4]=-2.;

  tgty[0]=0.3087;
  tgty[1]=0.3087;
  tgty[2]=.9087;
  tgty[3]=.9087;
  tgty[4]=.3087;
  auto target=new TPolyLine(5, tgtx, tgty);
  target->SetLineColor(kGreen);
  target->SetFillColor(kGreen);
  target->SetFillStyle(1001);
  //  target->SetTitle("target");
  target->Draw("f same");

auto scoperxgraph=new TGraph(3, scopex, scopey);

 auto scopeChanGraph=new TGraph();
  if(scopeChan<4&&scopeChan>=0){
    scopeChanGraph->SetPoint(0, scopex[scopeChan], scopey[scopeChan]);
    scopeChanGraph->SetMarkerSize(2.5);
    scopeChanGraph->SetMarkerStyle(4);
    scopeChanGraph->SetLineWidth(2);
    scopeChanGraph->SetMarkerColor(kViolet);
  }
 scoperxgraph->SetMarkerColor(kRed);
  scoperxgraph->SetMarkerStyle(20);
  scoperxgraph->SetMarkerSize(1.5);
  if(fScopeLoaded){
    scoperxgraph->Draw("p same");
    scopeChanGraph->Draw("p same");
  }
  
  auto surfrxgraph=new TGraph(12, surfx, surfy);
  surfrxgraph->SetMarkerColor(kRed);
  surfrxgraph->SetMarkerStyle(24);
  surfrxgraph->SetMarkerSize(1.5);
  if(fSurfLoaded){
    surfrxgraph->Draw("p same");
  }

  TLegend *leg = new TLegend(.83, .6, 1, .95);
  leg->AddEntry(txgraph, "TX", "p");
  leg->AddEntry(surfrxgraph, "SURF", "p");
  leg->AddEntry(scoperxgraph, "Scope", "p");
if(scopeChan<4&&scopeChan>=0){
    leg->AddEntry(scopeChanGraph, "CH"+TString::Itoa(scopeChan, 10), "p");
  }
  leg->AddEntry(target, "target", "l");
  leg->AddEntry(wall, "wall", "l");
  leg->AddEntry(support, "blocks", "l");
  leg->AddEntry(bp, "beam pipe", "l");
  
  leg->Draw();
  
  return 1;

  
}


TH2D* T576Event::pointingMap(double dx, int draw, int type){
  double tmin=20;
  double tmax=150;
  double dt[12][12]={0};
  double maxdt[12][12]={0};
  TGraph *grc[12][12]={new TGraph()};
  TVector3 source, d1, d2;
  vector<double>xx, yy,zz;
  auto graphs = vector<TGraph*>(12);
  auto surfx=vector<double>(12);
  auto surfy=vector<double>(12);
  double tot=0;
  auto deltat=1./(3.2*fInterpGSs);
  for(int i=0;i<12;i++){
    //graphs[i]=TUtil::normalize(TUtil::FFT::hilbertEnvelope(TUtil::getChunkOfGraph(surf->ch[i], tmin, tmax, 1)));

    if(type==2){
 graphs[i]=TUtil::normalize(TUtil::power((TUtil::getChunkOfGraph(surf->ch[i], tmin, tmax, 1))));
    }
    if(type==1){
      graphs[i]=TUtil::normalize(TUtil::FFT::hilbertEnvelope((TUtil::getChunkOfGraph(surf->ch[i], tmin, tmax, 1))));
    }
    else{
    graphs[i]=TUtil::normalize((TUtil::getChunkOfGraph(surf->ch[i], tmin, tmax, 1)));
    }
  }
  //auto window=rectangularWindow(10./deltat, 150./deltat, 250./deltat, deltat);
  for(int i=0;i<12;i++){
    for(int j=0;j<12;j++){
      double maxdelay=2.*(surf->pos[i]-surf->pos[j]).Mag()/TUtil::c_light;
      // cout<<i<<" "<<j<<" "<<maxdelay<<" "<<maxdelay/deltat<<endl;
      grc[i][j]=TUtil::crossCorrelate(graphs[i], graphs[j], maxdelay);
      // grc[i][j]=crossCorrelateUseful(graphs[j], graphs[i], graphs[j]->GetN(), (int)((graphs[i]->GetN()/2) -maxdt/deltat), (int)((graphs[i]->GetN()/2) +maxdt/deltat));
    }
  }
  //  setWarmPalette();
  double maxTot=0;
  TGraph *outgraphs[12];
  double coordincr=dx;
  TVector3 offset(.5, .5, .5);
  int nbins=16.01/dx;
  TH2D * gout=new TH2D("outhist", "outhist", nbins, -8, 8, nbins, -8, 8);
  double dtt = .2;
  for(double x=-8.01;x<8.01;x+=coordincr){

    source.SetZ((double)x);//z = x

    for(double y=-8.01;y<8.01;y+=coordincr){
      tot=0;
      source.SetX((double)y);//x=y

      xx.push_back((double)x);

      yy.push_back((double)y);
      for(int j=1;j<12;j++){
	//if(j<11&&j>2)continue;
	TVector3 surfupdated=surf->pos[j];//+offset;
	d1=source-surfupdated;//[j];
	//cout<<j<<" "<<i<<" "<<l<<endl;
	for(int k=j+1;k<12;k++){
	  //if(k<11&&k>2) continue;
	  //	  if(j==k)continue;
	  d2=source-surf->pos[k];
	  dt[j][k]=(abs(d2.Mag())-abs(d1.Mag()))/TUtil::c_light;
	  // cout<<dt[j][k]<<endl;
	  //	  tot+=grc[j][k]->Eval(dt[j][k]);
	  gout->Fill(x, y, grc[j][k]->Eval(dt[j][k])); 
	  //auto temp=TUtil::getChunkOfGraph(grc[j][k], dt[j][k]-dtt, dt[j][k]+dtt);
	  //tot+=TMath::MaxElement(temp->GetN(), temp->GetY());
	  //delete temp;

	}
      }

      zz.push_back(tot);
    }
  }


    //  auto gout=new TGraph2D(xx.size(), &xx[0], &yy[0], &zz[0]);
  gout->SetName("map"+minor);
  //    if(plot==1){
  //     auto crap=new TCanvas("map"+minor, "map"+minor, 700,600);

  //gout->SetMarkerStyle(20);
  gout->SetTitle("");

  gout->GetXaxis()->SetTitle("z (m)");
  gout->GetYaxis()->SetTitle("x (m)");
  //gout->GetHistogram()->GetXaxis()->SetTitle("z (m)");
  //gout->GetHistogram()->GetYaxis()->SetTitle("x (m)");

  if(draw==1){
    gout->SetMarkerStyle(20);
    gout->Draw("colz");
    //gPad->SetTheta(90.01);
    //gPad->SetPhi(0.01);
    //gPad->Update();


    //hist->Draw("colz");
    for(int i=0;i<12;i++){
      surfx[i]=surf->pos[i].Z();
      surfy[i]=surf->pos[i].X();
      //        cout<<surf[i].x<<" "<<surf[i].y<<" "<<surf[i].z<<endl;
    }
    auto rxgraph=new TGraph(12, &surfx[0], &surfy[0]);
    rxgraph->SetMarkerColor(kOrange);
    rxgraph->SetMarkerStyle(20);
    rxgraph->Draw("p same");

    //     auto rxtruegraph=new TGraph(12, surftruex, surftruey);
    //     rxtruegraph->SetMarkerColor(kRed);
    //     rxtruegraph->SetMarkerStyle(20);
    //     rxtruegraph->Draw("p same");

    auto txgraph=new TGraph();
    TVector3 tx(0,0,0);
    tx=txPos;
    txgraph->SetPoint(txgraph->GetN(), tx.Z(), tx.X());
    txgraph->SetMarkerColor(kYellow);
    txgraph->SetMarkerStyle(21);
    txgraph->Draw("p same");

    double tgtx[5], tgty[5];
    tgtx[0] = -2.;
    tgtx[1]=2.;
    tgtx[2]=2.;
    tgtx[3]=-2.;
    tgtx[4]=-2.;

    tgty[0]=0.3087;
    tgty[1]=0.3087;
    tgty[2]=.9087;
    tgty[3]=.9087;
    tgty[4]=.3087;
    auto target=new TPolyLine(5, tgtx, tgty);
    target->SetLineColor(kGreen);
    target->Draw("l same");

    double bpx[5], bpy[5];
    bpx[0] = -8.;
    bpx[1]=-3.6;
    bpx[2]=-3.6;
    bpx[3]=-8;

    bpy[0]=.5;
    bpy[1]=.5;
    bpy[2]=.7;
    bpy[3]=.7;
    bpy[4]=0;
    auto bp=new TPolyLine(4, bpx, bpy);
    bp->SetLineColor(kPink);
    bp->Draw("l same");

    double bdx[5], bdy[5];
    bdx[0] = 2;
    bdx[1]=5.6;
    bdx[2]=5.6;
    bdx[3]=2;
    bdx[4]=2;

    bdy[0]=2.44;
    bdy[1]=2.44;
    bdy[2]=-1.24;
    bdy[3]=-1.24;
    bdy[4]=2.44;
    auto bd=new TPolyLine(5, bdx, bdy);
    bd->SetLineColor(kPink);
    bd->Draw("l same");

    double wallx[5], wally[5];
    wallx[0] = 5.6;
    wallx[1]=5.6;

    wally[0]=8;
    wally[1]=-8;

    auto wall=new TPolyLine(2, wallx, wally);
    wall->SetLineColor(kPink);
    wall->Draw("l same");

    double bd2x[5], bd2y[5];
    bd2x[0] = 3.83;
    bd2x[1]=3.83;

    bd2y[0]=-1.24;
    bd2y[1]=2.44;

    auto bd2=new TPolyLine(2, bd2x, bd2y);
    bd2->SetLineColor(kPink);
    bd2->Draw("l same");

    double blockx[5], blocky[5];
    blockx[0] = -6;
    blockx[1]=-5;
    blockx[2]=-5.3;
    blockx[3]=-6.3;
    blockx[4]=-6;

    blocky[0]=-3;
    blocky[1]=-5;
    blocky[2]=-5.3;
    blocky[3]=-3.3;
    blocky[4]=-3.;
    auto block=new TPolyLine(5, blockx, blocky);
    block->SetLineColor(kPink);
    block->Draw("l same");

    double supportx[5], supporty[5];
    supportx[0] = -7.2075;
    supportx[1]=2;
    supportx[2]=2;
    supportx[3]=-7.2075;
    supportx[4]=-7.2075;

    supporty[0]=1.517;
    supporty[1]=1.517;
    supporty[2]=-.317;
    supporty[3]=-.317;
    supporty[4]=1.517;
    auto support=new TPolyLine(5, supportx, supporty);
    support->SetLineColor(kPink);
    support->Draw("l same");
  }
  return gout;
}


TH2D* T576Event::pointingMapDev(double dx, int draw, int hilbert, TVector3 * position){
  double tmin=20;
  double tmax=150;
  double dt[12][12]={0};
  double maxdt[12][12]={0};
  TGraph *grc[12][12]={new TGraph()};
  TVector3 source, d1, d2;
  vector<double>xx, yy,zz;
  auto graphs = vector<TGraph*>(12);
  auto surfx=vector<double>(12);
  auto surfy=vector<double>(12);
  double tot=0;
  auto deltat=1./(3.2*fInterpGSs);
  for(int i=0;i<12;i++){
    //graphs[i]=TUtil::normalize(TUtil::FFT::hilbertEnvelope(TUtil::getChunkOfGraph(surf->ch[i], tmin, tmax, 1)));
    
    if(hilbert==1){
      graphs[i]=TUtil::normalize(TUtil::FFT::hilbertEnvelope(TUtil::brickWallFilter(TUtil::getChunkOfGraph(surf->ch[i], tmin, tmax, 1), .9, 1.5)));
    }
    else{
    graphs[i]=TUtil::normalize(TUtil::brickWallFilter(TUtil::getChunkOfGraph(surf->ch[i], tmin, tmax, 1), .9, 1.5));
    }
  }
  //auto window=rectangularWindow(10./deltat, 150./deltat, 250./deltat, deltat);
  for(int i=0;i<12;i++){
    for(int j=0;j<12;j++){
      double maxdelay=(position[i]-position[j]).Mag()/TUtil::c_light;
      // cout<<i<<" "<<j<<" "<<maxdelay<<" "<<maxdelay/deltat<<endl;
      grc[i][j]=TUtil::crossCorrelate(graphs[i], graphs[j], maxdelay);
      // grc[i][j]=crossCorrelateUseful(graphs[j], graphs[i], graphs[j]->GetN(), (int)((graphs[i]->GetN()/2) -maxdt/deltat), (int)((graphs[i]->GetN()/2) +maxdt/deltat));
    }
  }
  //  setWarmPalette();
  double maxTot=0;
  TGraph *outgraphs[12];
  double coordincr=dx;
  TVector3 offset(.5, .5, .5);
  int nbins=16.01/dx;
  TH2D * gout=new TH2D("outhist", "outhist", nbins, -8, 8, nbins, -8, 8);
  double dtt = .2;
  for(double x=-8.01;x<8.01;x+=coordincr){

    source.SetZ((double)x);//z = x

    for(double y=-8.01;y<8.01;y+=coordincr){
      tot=0;
      source.SetX((double)y);//x=y

      xx.push_back((double)x);

      yy.push_back((double)y);
      for(int j=1;j<12;j++){
	//if(j<11&&j>2)continue;
	TVector3 surfupdated=position[j];//+offset;
	d1=source-surfupdated;//[j];
	//cout<<j<<" "<<i<<" "<<l<<endl;
	for(int k=j+1;k<12;k++){
	  //if(k<11&&k>2) continue;
	  //	  if(j==k)continue;
	  d2=source-position[k];
	  dt[j][k]=(abs(d1.Mag())-abs(d2.Mag()))/TUtil::c_light;
	  // cout<<dt[j][k]<<endl;
	  //	  tot+=grc[j][k]->Eval(dt[j][k]);
	  gout->Fill(x, y, grc[j][k]->Eval(dt[j][k])); 
	  //auto temp=TUtil::getChunkOfGraph(grc[j][k], dt[j][k]-dtt, dt[j][k]+dtt);
	  //tot+=TMath::MaxElement(temp->GetN(), temp->GetY());
	  //delete temp;

	}
      }

      zz.push_back(tot);
    }
  }


    //  auto gout=new TGraph2D(xx.size(), &xx[0], &yy[0], &zz[0]);
  gout->SetName("map"+minor);
  //    if(plot==1){
  //     auto crap=new TCanvas("map"+minor, "map"+minor, 700,600);

  //gout->SetMarkerStyle(20);
  gout->SetTitle("");

  gout->GetXaxis()->SetTitle("z (m)");
  gout->GetYaxis()->SetTitle("x (m)");
  //gout->GetHistogram()->GetXaxis()->SetTitle("z (m)");
  //gout->GetHistogram()->GetYaxis()->SetTitle("x (m)");

  if(draw==1){
    gout->SetMarkerStyle(20);
    gout->Draw("colz");
    //gPad->SetTheta(90.01);
    //gPad->SetPhi(0.01);
    //gPad->Update();


    //hist->Draw("colz");
    for(int i=0;i<12;i++){
      surfx[i]=position[i].Z();
      surfy[i]=position[i].X();
      //        cout<<surf[i].x<<" "<<surf[i].y<<" "<<surf[i].z<<endl;
    }
    auto rxgraph=new TGraph(12, &surfx[0], &surfy[0]);
    rxgraph->SetMarkerColor(kOrange);
    rxgraph->SetMarkerStyle(20);
    rxgraph->Draw("p same");

    //     auto rxtruegraph=new TGraph(12, surftruex, surftruey);
    //     rxtruegraph->SetMarkerColor(kRed);
    //     rxtruegraph->SetMarkerStyle(20);
    //     rxtruegraph->Draw("p same");

    auto txgraph=new TGraph();
    TVector3 tx(0,0,0);
    tx=txPos;
    txgraph->SetPoint(txgraph->GetN(), tx.Z(), tx.X());
    txgraph->SetMarkerColor(kYellow);
    txgraph->SetMarkerStyle(21);
    txgraph->Draw("p same");

    double tgtx[5], tgty[5];
    tgtx[0] = -2.;
    tgtx[1]=2.;
    tgtx[2]=2.;
    tgtx[3]=-2.;
    tgtx[4]=-2.;

    tgty[0]=0.3087;
    tgty[1]=0.3087;
    tgty[2]=.9087;
    tgty[3]=.9087;
    tgty[4]=.3087;
    auto target=new TPolyLine(5, tgtx, tgty);
    target->SetLineColor(kGreen);
    target->Draw("l same");

    double bpx[5], bpy[5];
    bpx[0] = -8.;
    bpx[1]=-3.6;
    bpx[2]=-3.6;
    bpx[3]=-8;

    bpy[0]=.5;
    bpy[1]=.5;
    bpy[2]=.7;
    bpy[3]=.7;
    bpy[4]=0;
    auto bp=new TPolyLine(4, bpx, bpy);
    bp->SetLineColor(kPink);
    bp->Draw("l same");

    double bdx[5], bdy[5];
    bdx[0] = 2;
    bdx[1]=5.6;
    bdx[2]=5.6;
    bdx[3]=2;
    bdx[4]=2;

    bdy[0]=2.44;
    bdy[1]=2.44;
    bdy[2]=-1.24;
    bdy[3]=-1.24;
    bdy[4]=2.44;
    auto bd=new TPolyLine(5, bdx, bdy);
    bd->SetLineColor(kPink);
    bd->Draw("l same");

    double wallx[5], wally[5];
    wallx[0] = 5.6;
    wallx[1]=5.6;

    wally[0]=8;
    wally[1]=-8;

    auto wall=new TPolyLine(2, wallx, wally);
    wall->SetLineColor(kPink);
    wall->Draw("l same");

    double bd2x[5], bd2y[5];
    bd2x[0] = 3.83;
    bd2x[1]=3.83;

    bd2y[0]=-1.24;
    bd2y[1]=2.44;

    auto bd2=new TPolyLine(2, bd2x, bd2y);
    bd2->SetLineColor(kPink);
    bd2->Draw("l same");

    double blockx[5], blocky[5];
    blockx[0] = -6;
    blockx[1]=-5;
    blockx[2]=-5.3;
    blockx[3]=-6.3;
    blockx[4]=-6;

    blocky[0]=-3;
    blocky[1]=-5;
    blocky[2]=-5.3;
    blocky[3]=-3.3;
    blocky[4]=-3.;
    auto block=new TPolyLine(5, blockx, blocky);
    block->SetLineColor(kPink);
    block->Draw("l same");

    double supportx[5], supporty[5];
    supportx[0] = -7.2075;
    supportx[1]=2;
    supportx[2]=2;
    supportx[3]=-7.2075;
    supportx[4]=-7.2075;

    supporty[0]=1.517;
    supporty[1]=1.517;
    supporty[2]=-.317;
    supporty[3]=-.317;
    supporty[4]=1.517;
    auto support=new TPolyLine(5, supportx, supporty);
    support->SetLineColor(kPink);
    support->Draw("l same");
  }
  return gout;
}




TVector3* T576Event::fixPositionsDev(double dx, int maxIter, int hilbert,TVector3 source, TVector3 *positions){
  double tmin=20;
  double tmax=150;
  double dt[12][12]={0};
  double maxdt[12][12]={0};
  TGraph *grc[12][12]={new TGraph()};
  TVector3 d1, d2;
  vector<double>xx, yy,zz;
  auto graphs = vector<TGraph*>(12);
  auto surfx=vector<double>(12);
  auto surfy=vector<double>(12);

  
  // for(int i=0;i<12;i++){
  //   positions=(TVector3)surf->pos[i].Clone();
  // }
  double tot=0;
  auto deltat=1./(3.2*fInterpGSs);
  for(int i=0;i<12;i++){
    //graphs[i]=TUtil::normalize(TUtil::FFT::hilbertEnvelope(TUtil::getChunkOfGraph(surf->ch[i], tmin, tmax, 1)));
    
    if(hilbert==1){
      graphs[i]=TUtil::normalize(TUtil::FFT::hilbertEnvelope(TUtil::bandpassFilter(TUtil::getChunkOfGraph(surf->ch[i], tmin, tmax, 1), .9, 1.5)));
    }
    else{
    graphs[i]=TUtil::normalize(TUtil::bandpassFilter(TUtil::getChunkOfGraph(surf->ch[i], tmin, tmax, 1), .9, 1.5));
    }
  }
  //auto window=rectangularWindow(10./deltat, 150./deltat, 250./deltat, deltat);
  for(int i=0;i<12;i++){
    for(int j=0;j<12;j++){
      //to give a larger window of coincidence to move the postions
      double maxdelay=(positions[i]-positions[j]).Mag()+5./TUtil::c_light;
      // cout<<i<<" "<<j<<" "<<maxdelay<<" "<<maxdelay/deltat<<endl;
      grc[i][j]=TUtil::crossCorrelate(graphs[i], graphs[j], maxdelay);
      // grc[i][j]=crossCorrelateUseful(graphs[j], graphs[i], graphs[j]->GetN(), (int)((graphs[i]->GetN()/2) -maxdt/deltat), (int)((graphs[i]->GetN()/2) +maxdt/deltat));
    }
  }
  //  setWarmPalette();
  double maxTot=0;
  TGraph *outgraphs[12];
  double coordincr=dx;
  TVector3 offset(.5, .5, .5);
  int nbins=16.01/dx;
  TH2D * gout=new TH2D("outhist", "outhist", nbins, -8, 8, nbins, -8, 8);
  double dtt = .2;

  // double x = .6;
  // double y = -2.;

  // source.SetZ((double)x);//z = x
  // source.SetX((double)y);//x=y
  double lastTot=0;
  int nIter=0;
  //  double ddx=.01;
  TVector3 initPos;
  while(nIter<maxIter){

    for(int j=1;j<12;j++){

      for(int i=0;i<4;i++){
	initPos=positions[j];
	switch(i){
	case(0):
	  positions[j].SetMag(positions[j].Mag()+dx);
	  break;
	case(1):
	  positions[j].SetMag(positions[j].Mag()-dx);
	  break;
	case(2):
	  positions[j].SetTheta(positions[j].Theta()+(dx));
	  break;
	case(3):
	  positions[j].SetTheta(positions[j].Theta()-(dx));
	  break;
	}
      
	//if(j<11&&j>2)continue;
	d1=source-positions[j];//[j];
	//cout<<j<<" "<<i<<" "<<l<<endl;
	for(int k=1;k<12;k++){
	  //if(k<11&&k>2) continue;
	  //	  if(j==k)continue;
	  d2=source-positions[k];
	  dt[j][k]=(abs(d1.Mag())-abs(d2.Mag()))/TUtil::c_light;
	  // cout<<dt[j][k]<<endl;
	  //	  tot+=grc[j][k]->Eval(dt[j][k]);
	  tot+=grc[j][k]->Eval(dt[j][k]);
	  //	  gout->Fill(x, y, grc[j][k]->Eval(dt[j][k])); 
	  //auto temp=TUtil::getChunkOfGraph(grc[j][k], dt[j][k]-dtt, dt[j][k]+dtt);
	  //tot+=TMath::MaxElement(temp->GetN(), temp->GetY());
	  //delete temp;

	}
	positions[j]=tot>lastTot?positions[j]:initPos;

	lastTot=tot;
	tot=0.;
      }
      
    }
    //    cout<<positions[1].X()<<endl;
    nIter++;
  }
  
  


  return positions;
}


vector<TPolyLine*> T576Event::getTheRoom(Color_t lineColor, Color_t tgtColor){
    double tgtx[5], tgty[5];
    tgtx[0] = -2.;
    tgtx[1]=2.;
    tgtx[2]=2.;
    tgtx[3]=-2.;
    tgtx[4]=-2.;

    tgty[0]=0.3087;
    tgty[1]=0.3087;
    tgty[2]=.9087;
    tgty[3]=.9087;
    tgty[4]=.3087;
    auto target=new TPolyLine(5, tgtx, tgty);
    target->SetLineColor(tgtColor);
    target->Draw("l same");

    double bpx[5], bpy[5];
    bpx[0] = -8.;
    bpx[1]=-3.6;
    bpx[2]=-3.6;
    bpx[3]=-8;

    bpy[0]=.5;
    bpy[1]=.5;
    bpy[2]=.7;
    bpy[3]=.7;
    bpy[4]=0;
    auto bp=new TPolyLine(4, bpx, bpy);
    bp->SetLineColor(lineColor);
    bp->Draw("l same");

    double bdx[5], bdy[5];
    bdx[0] = 2;
    bdx[1]=5.6;
    bdx[2]=5.6;
    bdx[3]=2;
    bdx[4]=2;

    bdy[0]=2.44;
    bdy[1]=2.44;
    bdy[2]=-1.24;
    bdy[3]=-1.24;
    bdy[4]=2.44;
    auto bd=new TPolyLine(5, bdx, bdy);
    bd->SetLineColor(lineColor);
    bd->Draw("l same");

    double wallx[5], wally[5];
    wallx[0] = 5.6;
    wallx[1]=5.6;

    wally[0]=8;
    wally[1]=-8;

    auto wall=new TPolyLine(2, wallx, wally);
    wall->SetLineColor(lineColor);
    wall->Draw("l same");

    double bd2x[5], bd2y[5];
    bd2x[0] = 3.83;
    bd2x[1]=3.83;

    bd2y[0]=-1.24;
    bd2y[1]=2.44;

    auto bd2=new TPolyLine(2, bd2x, bd2y);
    bd2->SetLineColor(lineColor);
    bd2->Draw("l same");

    double blockx[5], blocky[5];
    blockx[0] = -6;
    blockx[1]=-5;
    blockx[2]=-5.3;
    blockx[3]=-6.3;
    blockx[4]=-6;

    blocky[0]=-3;
    blocky[1]=-5;
    blocky[2]=-5.3;
    blocky[3]=-3.3;
    blocky[4]=-3.;
    auto block=new TPolyLine(5, blockx, blocky);
    block->SetLineColor(lineColor);
    block->Draw("l same");

    double supportx[5], supporty[5];
    supportx[0] = -7.2075;
    supportx[1]=2;
    supportx[2]=2;
    supportx[3]=-7.2075;
    supportx[4]=-7.2075;

    supporty[0]=1.517;
    supporty[1]=1.517;
    supporty[2]=-.317;
    supporty[3]=-.317;
    supporty[4]=1.517;
    auto support=new TPolyLine(5, supportx, supporty);
    support->SetLineColor(lineColor);
    support->Draw("l same");



    vector<TPolyLine*> lines;
    
    lines.push_back(target);
    lines.push_back(wall);
    lines.push_back(bd);
    lines.push_back(bd2);
    lines.push_back(block);
    lines.push_back(support);
    return lines;
}

TGraph * T576Event::getAntennas(Color_t color){
     //hist->Draw("colz");
  auto surfx=vector<double>(12);
  auto surfy=vector<double>(12);
  for(int i=0;i<12;i++){
      surfx[i]=surf->pos[i].Z();
      surfy[i]=surf->pos[i].X();
      //        cout<<surf[i].x<<" "<<surf[i].y<<" "<<surf[i].z<<endl;
    }
    auto rxgraph=new TGraph(12, &surfx[0], &surfy[0]);
    rxgraph->SetMarkerColor(color);
    rxgraph->SetMarkerStyle(20);
    rxgraph->Draw("p same");

  return rxgraph;
}

double T576Event::getInterpGSs(){
  return fInterpGSs;
}


TNtuple * T576Event::integrateAllWithSideband(int major, int minor, int scopeOrSurf, int channel, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, double sbxmin, double sbxmax, double sbymin, double sbymax, double scale){
  int number=0;
  TNtuple *tup=new TNtuple("tup", "tup", "sig:sb");
  auto graphs=vector<TGraph*>();
  if(scopeOrSurf==0){
    //    loadScopeEvent(major, minor, 0);
    for(int i=0;i<scopeNEvents;i++){
      loadScopeEvent(major, minor, i);
      if(isGood){
	auto spec = TUtil::FFT::spectrogram(scope->ch[channel], nfft, overlap, zeroPadLength, window, dbFlag);
	spec->Scale(scale);
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

TNtuple * T576Event::sidebandSubtractAll(int major, int minor, int scopeOrSurf, int channel, int nfft, int overlap, int zeroPadLength, int window, int dbFlag, double xmin, double xmax, double ymin, double ymax, double scale){
  int number=0;
  TNtuple *tup=new TNtuple("tup", "tup", "sig:err");
  auto graphs=vector<TGraph*>();
  if(scopeOrSurf==0){
    //    loadScopeEvent(major, minor, 0);
    for(int i=0;i<scopeNEvents;i++){
      loadScopeEvent(major, minor, i);
      if(isGood){
	auto spec = TUtil::FFT::spectrogram(scope->ch[channel], nfft, overlap, zeroPadLength, window, dbFlag);
	double err=0.;
	spec->Scale(scale);
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

