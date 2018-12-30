#include "t576/T576Event.hh"

int main(){
  time_t time_start, time_end;

  auto ev=new T576Event();
  //  ev->setInterpGsNs(20.);
  int eventCount=0;
  cout<<"running through first several thousand events. skipping run 0 because it was a test run..."<<endl;
  for(int i=0;i<10000;i++){
    if(i==1)time(&time_start);
  
    eventCount+=ev->loadScopeEvent(i);
    cout.flush()<<i<<"     \r";
  
  }

  time(&time_end);
  cout.flush();

  printf("read first %i events in %li seconds.\n", eventCount, time_end-time_start); 
  delete(ev);
  exit(0);
}
  
