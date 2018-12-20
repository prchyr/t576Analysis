#include "t576/T576Event.hh"

int main(){
  time_t time_start, time_end;

  auto ev=new T576Event();
  int eventCount=0;
  for(int i=0;i<50000;i++){
    if(i==1)time(&time_start);
    eventCount+=ev->loadScopeEvent(i);
    cout.flush()<<i<<"     \r";
  }
  ev->scope->gr[3]->Draw("al");
  time(&time_end);
  cout.flush();
  //  cout<<time_end<<" "<<time_start<<endl;
  printf("read first %i events in %li seconds.\n", eventCount, time_end-time_start); 
  delete(ev);
  exit(0);
}
  
