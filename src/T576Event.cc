#include "T576Event.hh"

using namespace CLHEP;
using namespace std;

// T576Event::T576Event(int run_major, int run_minor,int event){
//   loadEvent(int run_major, int run_minor,int event);
// }


int T576Event::loadEvent(int run_major, int run_minor,int event){
  TString major=TString::Itoa(run_major, 10);
  TString minor=TString::Itoa(run_minor, 10);
  
  TString top_dir = getenv("T576_DATA_DIR");
  TString directory=top_dir+"/root/";
  
  if(!directory){
    cout<<"T576_DATA_DIR not set. please set this flag so that the data can be found. this should be the top directory inside of which is py/ and root/."<<endl;
    return (-1);
  }
  
  TString ends_with=major+"_"+minor;
  TString filename=fFilename;
  

  //check that this file isn't already loaded, if it is, save some time
  if(filename.EndsWith(ends_with+".root")){
    goto skip;
  }
  TSystemDirectory dir(directory, directory);
  TList *files = dir.GetListOfFiles();

  //find the filename
  if(files){
    TSystemFile *file;
    TString fname;
    TString ext = ".root";
    TIter next(files);
    while((file=(TSystemFile*)next())){
      fname = file->GetName();
      if(!file->IsDirectory()&&fname.EndsWith(ends_with+ext)){
	filename=fname;
	//	cout<<filename<<endl;
      }
    }
    delete (file);
  }
  delete (files);
  if(!filename){
    cout<<"file not found!"<<endl;
    return (-1);
  }


 skip:
  auto file=TFile::Open(directory+filename);
  auto tree=(TTree*)file->Get("tree");
  

  tree->SetBranchAddress("ch1", scope.ch[0]);
  tree->SetBranchAddress("ch2", scope.ch[1]);
  tree->SetBranchAddress("ch3", scope.ch[2]);
  tree->SetBranchAddress("ch4", scope.ch[3]);
  tree->SetBranchAddress("time", scope.time);
  tree->SetBranchAddress("timestamp", &timestamp);
  tree->GetEntry(event);
  
  auto length=sizeof(scope.time)/sizeof(*scope.time);

  scope.time[0]=0.;
  scope.time[19999]=scope.time[19998]+.05;
  for(int i=0;i<4;i++){
    auto graph=new TGraph(length, scope.time, scope.ch[i]);
    *scope.gr[i]=*graph;
    delete(graph);
  }

  delete(tree);

  return 1;
}

//constructor with a tree

int T576Event::loadEvent(int run_major, int run_minor,TTree * tree){
  //double * ch[4];
  //  double * time;
  //ULong64_t timestamp;
  tree->SetBranchAddress("ch1", scope.ch[0]);
  tree->SetBranchAddress("ch2", scope.ch[1]);
  tree->SetBranchAddress("ch3", scope.ch[2]);
  tree->SetBranchAddress("ch4", scope.ch[3]);
  tree->SetBranchAddress("time", scope.time);
  tree->SetBranchAddress("timestamp", &timestamp);

  auto length=sizeof(scope.time)/sizeof(*scope.time);

  for(int i=0;i<4;i++){
    auto graph=new TGraph(length, scope.time, scope.ch[i]);
    *scope.gr[i]=*graph;
    delete(graph);
  }

  //  charge=getCharge(scope.ch[3]);
  
}
