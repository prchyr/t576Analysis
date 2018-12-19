#include "t576/T576Event.hh"

int main(){
  auto ev=new T576Event();
  for(int i=0;i<10000;i++){
    ev->loadScopeEvent(i);
    cout.flush()<<i<<"     \r";
  }
  delete(ev);
  exit(0);
}
  
