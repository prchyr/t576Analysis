/*
copyright s. prohira and the T576 collaboration 2018
released under the GNU General Public License version 3
*/

#ifndef TUTILGRAPH_H
#define TUTILGRAPH_H

#include "TGraph.h"

//using namespace CLHEP;
//using namespace std;
/*
This class is just so we can use some operators and things to make life simpler.

re-implements the tgraph constructors. 
 */

class TUtilGraph: public TGraph{
  public:
  TUtilGraph(Int_t n);
  TUtilGraph(Int_t n, const Double_t * x, const Double_t * y);
  TUtilGraph(const TGraph *gr);
  
  double sum();
  double sumPower();
  double rms();
  double mean();
  
  virtual ~TUtilGraph(){};

  TUtilGraph * operator*(const TUtilGraph *b);
  TUtilGraph * operator+(const TUtilGraph *b);
  TUtilGraph * operator-(const TUtilGraph *b);
  TUtilGraph * operator*(const double b);
  TUtilGraph * operator+(const double b);
  TUtilGraph * operator-(const double b);
  
protected:
    ClassDef(TUtilGraph, 1);
};

#endif
