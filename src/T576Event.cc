#include "T576Event.hh"

using namespace CLHEP;
using namespace std;




//on construction, check that the data directory structure is good and indexed. if not, build an index file for quick grabbing of events.
int T576Event::checkStatus(){
  fInstallDir=getenv("T576_INSTALL_DIR");
  if(fInstallDir==""){
    fInstallDir="/usr/local";
  }

  TFile* indexFile;
  indexFile=TFile::Open(fInstallDir+"/share/t576/eventIndex.root");

  if(!indexFile){
    cout<<"event index not built yet. building..."<<endl;
    buildEventIndex();
    indexFile=TFile::Open(fInstallDir+"/share/t576/eventIndex.root");
  }
  
  fIndexTree=(TTree*)indexFile->Get("indexTree");
  fIndexTree->SetBranchAddress("scopeFilename", &scopeFilename);
  fIndexTree->SetBranchAddress("major", &major);
  fIndexTree->SetBranchAddress("minor", &minor);
  fIndexTree->SetBranchAddress("subEvNo", &subEvNo);
  fIndexTree->SetBranchAddress("scopeEvNo", &scopeEvNo);
  fIndexTree->SetBranchAddress("surfEvNo", &surfEvNo);
  
  fIndexBuilt=1;
  fIndexTree->BuildIndex("scopeEvNo", "scopeEvNo");
  fScopeEvNoIndex=(TTreeIndex*)fIndexTree->GetTreeIndex();
  fIndexTree->BuildIndex("major", "minor");
  fMajorMinorIndex=(TTreeIndex*)fIndexTree->GetTreeIndex();

}

//load a scope event using a run major/minor combination and event within that file

int T576Event::loadScopeEvent(int run_major, int run_minor,int event){

  //TString major=TString::Itoa(run_major, 10);
  //TString minor=TString::Itoa(run_minor, 10);
  
  if(fIndexBuilt==0){
    cout<<"event index not built yet. building..."<<endl;
    buildEventIndex();
  }
  
  // TString install_dir=getenv("T576_INSTALL_DIR");
  // if(install_dir==""){
  //   install_dir="/usr/local";
  // }

  // TFile* indexFile;
  // indexFile=TFile::Open(install_dir+"/share/t576/eventIndex.root");

  // if(indexFile->IsZombie()){
  //   cout<<"event index not built yet. building..."<<endl;
  //   buildEventIndex();
  //   indexFile=TFile::Open(install_dir+"/share/t576/eventIndex.root");
  // }

  TString top_dir = getenv("T576_DATA_DIR");
  TString directory=top_dir+"/root/";
  if(!directory){
    cout<<"T576_DATA_DIR not set. please set this flag so that the data can be found. this should be the top directory inside of which is py/ and root/."<<endl;
    return (-1);
  }
  
  //cout<<"asdf"<<endl;
  //  fIndexTree->SetBranchAddress("scopeFilename", &scopeFilename);
  fIndexTree->SetTreeIndex(fMajorMinorIndex);
  fIndexTree->GetEntry(fIndexTree->GetEntryNumberWithBestIndex(run_major, run_minor));
  //  fIndexTree->GetEntry(0);
  //  cout<<scopeFilename->Data();
  TString thisScopeFilename=scopeFilename->Data();
  if(thisScopeFilename==fScopeFilename)goto fileSkip;

  fScopeFilename=thisScopeFilename;
  //open the file
  //  cout<<scopeScopeFilename<<endl;
  fEventFile=TFile::Open(directory+thisScopeFilename);
  fEventTree=(TTree*)fEventFile->Get("tree");
  

  fEventTree->SetBranchAddress("ch1", scope->ch[0]);
  fEventTree->SetBranchAddress("ch2", scope->ch[1]);
  fEventTree->SetBranchAddress("ch3", scope->ch[2]);
  fEventTree->SetBranchAddress("ch4", scope->ch[3]);
  fEventTree->SetBranchAddress("time", scope->time);
  fEventTree->SetBranchAddress("timestamp", &timestamp);


 fileSkip:
  if(event>fEventTree->GetEntries()){
    cout<<"event number too high! this major/minor combination only contains "<<fEventTree->GetEntries()<<" events."<<endl<<"select another event number, or call loadScopeEvent(event) or loadSurfEvent(event) without major/minor to use overall event index"<<endl;
    return (-1);
  }

  fEventTree->GetEntry(event);

  //check the length of the record.
  auto length=sizeof(scope->time)/sizeof(*scope->time);

  //fill the event graphs for the scope->
  //fix the first and last values, which were recorded incorrectly
  scope->time[0]=0.;
  scope->time[19999]=scope->time[19998]+.05;
  for(int i=0;i<4;i++){
    auto graph=new TGraph(length, scope->time, scope->ch[i]);
    *scope->gr[i]=*graph;
    delete(graph);
  }
  
  //delete(tree);
  fEventFile->Close();
  //delete(file);
  getCharge(scope->gr[3]);
  
  //  delete (files);
  
  
  return 1;
  }


//load an event from the global scope event number index 
int T576Event::loadScopeEvent(int event){

  if(fIndexBuilt==0){
    cout<<"event index not built yet. building..."<<endl;
    buildEventIndex();
  }

  TString top_dir = getenv("T576_DATA_DIR");
  TString directory=top_dir+"/root/";
  if(!directory){
    cout<<"T576_DATA_DIR not set. please set this flag so that the data can be found. this should be the top directory inside of which is py/ and root/."<<endl;
    return (-1);
  }
  
  //cout<<"asdf"<<endl;
  //  fIndexTree->SetBranchAddress("scopeFilename", &scopeFilename);
  fIndexTree->SetTreeIndex(fScopeEvNoIndex);
  fIndexTree->GetEntry(fIndexTree->GetEntryNumberWithBestIndex(event, event));
  //  fIndexTree->GetEntry(0);
  //  cout<<thisScopeFilename->Data();
  TString thisScopeFilename=scopeFilename->Data();
  //open the file
  //  cout<<thisScopeFilename<<endl<<" "<<fScopeFilename<<endl;
  if(thisScopeFilename==fScopeFilename)goto fileSkip;

  cout<<thisScopeFilename<<endl;
  fEventFile=TFile::Open(directory+thisScopeFilename);
  fEventTree=(TTree*)fEventFile->Get("tree");

  //  memset(scope->ch, 0, sizeof(scope->ch));
  //memset(scope->gr, 0, sizeof(scope->gr));

  fEventTree->SetBranchAddress("ch1", scope->ch[0]);
  fEventTree->SetBranchAddress("ch2", scope->ch[1]);
  fEventTree->SetBranchAddress("ch3", scope->ch[2]);
  fEventTree->SetBranchAddress("ch4", scope->ch[3]);
  fEventTree->SetBranchAddress("time", scope->time);
  fEventTree->SetBranchAddress("timestamp", &timestamp);

  fScopeFilename=thisScopeFilename;
  fileSkip:
  if(subEvNo>=fEventTree->GetEntries()){
    cout<<"event number too high! this major/minor combination only contains "<<fEventTree->GetEntries()<<" events."<<endl<<"select another event number, or call loadEvent(event) without major/minor to use overall event index"<<endl;
  }

  fEventTree->GetEntry(event);

  //check the length of the record.
  auto length=sizeof(scope->time)/sizeof(*scope->time);

  //fill the event graphs for the scope->
  //fix the first and last values, which were recorded incorrectly
  scope->time[0]=0.;
  scope->time[19999]=scope->time[19998]+.05;
  for(int i=0;i<4;i++){
    auto graph=new TGraph(length, scope->time, scope->ch[i]);
    *scope->gr[i]=*graph;
    delete(graph);
  }
  
  //delete(tree);

  if(subEvNo==fEventTree->GetEntries())fEventFile->Close();
  //delete(file);
  getCharge(scope->gr[3]);
  
  //  delete (files);
  
  
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

int T576Event::Scope::getAntennaPositions(int run_major, int run_minor){

  
}

int T576Event::Surf::getAntennaPositions(int run_major, int run_minor){


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
  TString scopeFilename;
  int maj=0, min=0, subEvNo=0, scopeEvNo=0, surfEvNo=0;
  ULong64_t tstamp;
  indexTree->Branch("scopeFilename", &scopeFilename);
  indexTree->Branch("major", &maj);
  indexTree->Branch("minor", &min);
  indexTree->Branch("subEvNo", &subEvNo);
  indexTree->Branch("scopeEvNo", &scopeEvNo);
  indexTree->Branch("surfEvNo", &surfEvNo);
  indexTree->Branch("tstamp", &tstamp);

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
	      tree->SetBranchAddress("timestamp", &tstamp);
	      int nentries=tree->GetEntries();
	      cout.flush()<<fname<<"                  \r";	      
	      for(int j=0;j<nentries;j++){
		//cout<<scopeFilename<<" "<<major<<" "<<minor<<" "<<subEvNo<<" "<<scopeEvNo<<endl;
		tree->GetEntry(j);
		subEvNo=j;
		//	      scopeFilename=const_cast<char*>(fname.Data());
		scopeFilename=fname;
		indexTree->Fill();
		scopeEvNo++;
	      }
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
    return -1;
  }
  index->Write();
  index->Close();
  delete(index);
  delete(indexTree);

  fIndexBuilt=1;
  
  return 1;
}
