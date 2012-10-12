#include <iostream>


#ifndef PROGRESS_BAR
#define PROGRESS_BAR

struct PBar
{
  PBar(int i)   : progress(i){}
  PBar(double i): progress(i*100.0){}
  
  friend ostream& operator << (ostream& out, const PBar& pb ){
    out << pb.progress << flush;
    return out;
  }
  
  double progress;
};

#endif

