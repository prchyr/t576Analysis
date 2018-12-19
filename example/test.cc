#include "t576/T576Event.hh"

int main(){
  time_t time_start, time_end;
  time(&time_start);
  auto ev=new T576Event();
  for(int i=0;i<10000;i++){
    ev->loadScopeEvent(i);
    cout.flush()<<i<<"     \r";
  }
  time(&time_end);
  cout.flush();
  //  cout<<time_end<<" "<<time_start<<endl;
  printf("read first 10k events in %d seconds.\n", time_end-time_start); 
  delete(ev);
  exit(0);
}
  
