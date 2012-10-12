#include <iostream>


#ifndef FORMAT_OUTPUT
#define FORMAT_OUTPUT

struct PBar
{
  PBar(int i)   : progress(i){}
  PBar(double i): progress(i*100.0){}
  
  friend ostream& operator << (ostream& out, const PBar& pb ){
    for (int i=0;i<100;i++)
      out << '\b';
    out << pb.progress << '%' << flush;
    return out;
  }
  
  double progress;
};

#endif

