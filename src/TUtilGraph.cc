/*
copyright s. prohira and the T576 collaboration 2018
released under the GNU General Public License version 3
*/

#include "TUtilGraph.hh"
#include "TROOT.h"

ClassImp(TUtilGraph);

TUtilGraph::TUtilGraph(Int_t n)
{
  fNpoints = n;
  if (!CtorAllocate()) return;
  FillZero(0, fNpoints);
}

TUtilGraph::TUtilGraph(Int_t n, const Double_t *x, const Double_t *y)
{
  if (!x || !y) {
    fNpoints = 0;
  } else {
    fNpoints = n;
  }
  if (!CtorAllocate()) return;
  n = fNpoints * sizeof(Double_t);
  memcpy(fX, x, n);
  memcpy(fY, y, n);
}

TUtilGraph::TUtilGraph(const TGraph *gr)
{
  double * x = gr->GetX();
  double * y = gr->GetY();
  int n = gr->GetN();
  if (!x || !y) {
    fNpoints = 0;
  } else {
    fNpoints = n;
  }
  if (!CtorAllocate()) return;
  n = fNpoints * sizeof(Double_t);
  memcpy(fX, x, n);
  memcpy(fY, y, n);
}

double TUtilGraph::sum(){
  double val=0;
  for(int i=0;i<this->GetN();i++)val+=this->GetY()[i];
  return val;
}

double TUtilGraph::sumPower(){
  double val=0;
  for(int i=0;i<this->GetN();i++)val+=(this->GetY()[i]*this->GetY()[i]);
  return val;
}

double TUtilGraph::mean(){
  double val=0;
  for(int i=0;i<this->GetN();i++)val+=this->GetY()[i];
  return val/this->GetN();
}

double TUtilGraph::rms(){
  double val=0;
  for(int i=0;i<this->GetN();i++)val+=(this->GetY()[i]*this->GetY()[i]);
  return sqrt(val/this->GetN());
}


TUtilGraph * TUtilGraph::operator*(const TUtilGraph *b){
  for(int i=0;i<this->GetN();i++)this->GetY()[i]*=b->GetY()[i];
}


TUtilGraph * TUtilGraph::operator+(const TUtilGraph *b){
  for(int i=0;i<this->GetN();i++)this->GetY()[i]+=b->GetY()[i];

}


TUtilGraph * TUtilGraph::operator-(const TUtilGraph *b){
  auto thing=(TUtilGraph*)this->Clone();
  for(int i=0;i<thing->GetN();i++){
    thing->GetY()[i]=this->GetY()[i]-b->GetY()[i];
  }
  return thing;
}


TUtilGraph * TUtilGraph::operator*(const double b){
  for(int i=0;i<this->GetN();i++)this->GetY()[i]*=b;
}


TUtilGraph * TUtilGraph::operator+(const double b){
  for(int i=0;i<this->GetN();i++)this->GetY()[i]+=b;
}


TUtilGraph * TUtilGraph::operator-(const double b){
  for(int i=0;i<this->GetN();i++)this->GetY()[i]-=b;
}
