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
}











//load a scope event using a run major/minor combination and event within that file

int T576Event::loadScopeEvent(int run_major, int run_minor,int event){
  if(fIndexBuilt==0){
    cout<<"event index not built yet. building..."<<endl;
    buildEventIndex();
  }
  TString top_dir = getenv("T576_DATA_DIR");

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
  if(loadScopeEvent(scopeEvNo+event)==1){
    getCharge(scope->ch[3]);
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
int T576Event::loadScopeEvent(int event){

  if(fIndexBuilt==0){
    cout<<"event index not built yet. building..."<<endl;
    buildEventIndex();
  }
  
  TString top_dir = getenv("T576_DATA_DIR");

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
  //cout<<thisScopeFilename<<endl<<" "<<fScopeFilename<<endl;
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
  auto length=sizeof(scope->time)/sizeof(*scope->time);
  if(length!=20000)cout<<length;
  //cout<<major<<" "<<minor<<" "<<subEvNo<<" "<<length<<" "<<fEventTree->GetEntries()<<" "<<scope->dat[1][18]<<endl;
  //fill the event graphs for the scope->
  //fix the first and last values, which were recorded incorrectly
  scope->time[0]=0.;
  scope->time[19999]=scope->time[19998]+.05;
  for(int i=0;i<4;i++){
    TGraph * graph=new TGraph(length, scope->time, scope->dat[i]);
    graph->SetTitle("");
    graph->SetName("ch"+TString::Itoa(i, 10));
    graph->GetXaxis()->SetTitle("Time (ns)");
    graph->GetYaxis()->SetTitle("Volts (V)");
    graph->GetYaxis()->SetTitleOffset(1.15);
    graph->GetHistogram()->SetName("");
    if(fInterpGSs>0.){
      getInterpolatedGraph(graph, scope->ch[i]);
    }
    else{
      *scope->ch[i]=*graph;
    }
   delete(graph);
  }
  //  cout<<"here"<<endl;
  if(subEvNo==fEventTree->GetEntries())fEventFile->Close();
  //delete(file);
  getCharge(scope->ch[3]);
  
  //  delete (files);
  
  
  return 1;
}



int T576Event::loadSurfEvent(int run_major, int run_minor, int event){

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
  if(surfFilename->Data()=="null"){
    cout<<"no surf events for this run."<<endl;
    return 0;
  }
  if(major<1)return 0;
  if(event>surfNEvents){
    cout<<"event number too large for this major/minor combination"<<endl;
    return 0;
  }

  if(loadSurfEvent(surfEvNo+event)==1)return 1;
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
int T576Event::loadSurfEvent(int event){

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
    surf->pos[i].SetMagThetaPhi(surf->dist[i], surf->ang[i]*pi/180., pi);
  }
  if(major<1)return 0;

  TString thisSurfFilename=surfFilename->Data();
  //open the file
  //cout<<thisSurfFilename<<endl<<" "<<event<<endl;
  if(thisSurfFilename!=fSurfFilename){
    //load the event file
    TString openf=directory+thisSurfFilename;

    fDataset = cnpy::npz_load(openf.Data());

    //load the data. it is stored as an array of shorts of length (N events) x (12 channels ) x (1024 samples)

    fSurfDataArray = fDataset["data"];

    fSurfData=fSurfDataArray.data<short>();

    cnpy::NpyArray times_arr = fDataset["times"];
    fSurfTimes = times_arr.data<double>();

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
    surf->delays[i]=(surf->cableLengths[i]/3.2808)/(c_light*surf->velocityFactor);

    TGraph * graph=new TGraph(len, TUtil::makeIndices(len, 1./3.2, surf->delays[i]), surf->dat[i]);

    TGraph *grChunk=TUtil::getChunkOfGraph(graph, 0, 250);
    grChunk->SetTitle("");
    grChunk->SetName("ch"+TString::Itoa(i, 10));
    grChunk->GetXaxis()->SetTitle("Time (ns)");
    grChunk->GetYaxis()->SetTitle("Volts (V)");
    grChunk->GetYaxis()->SetTitleOffset(1.15);
    grChunk->GetHistogram()->SetName("");
    if(fInterpGSs>0.){
      getInterpolatedGraph(grChunk, surf->ch[i]);
    }
    else{
      *surf->ch[i]=*grChunk;
    }
    delete(grChunk);
    delete(graph);
  }
  
  //getCharge();  
  return 1;
}


double T576Event::getCharge(TGraph *ict){
  double tot=0.;
  for(int i=0;i<ict->GetN();i++){
    tot+=ict->GetY()[i];
  }
  double xval=ict->GetX()[10]-ict->GetX()[9];

  charge = .4*tot*xval;
  return charge;
}


int T576Event::getInterpolatedGraph(TGraph * inGraph, TGraph *outGraph){
  ROOT::Math::Interpolator interp(inGraph->GetN(), ROOT::Math::Interpolation::kAKIMA);
  interp.SetData(inGraph->GetN(), inGraph->GetX(), inGraph->GetY());

  //get dt, assuming even sampling.
  double inDt=inGraph->GetX()[50]-inGraph->GetX()[49];
  double inGsPerS=1./inDt;

  double outDt=1./fInterpGSs;

  int samps=(int) (inGraph->GetN()*(fInterpGSs/inGsPerS));

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
