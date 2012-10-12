#include "RanlibCPP.h"
#include <iostream>

using namespace std;


int main()
{
//    char  alphabet[] = {'a','b','c','d','e'};

  
//   for (int i=0;i<10;i++){ // make 10 sequences
//     cout << "$\\mathbf{";
//     for (int j=0;j<100;j++){ // of length 100
//       char l = alphabet[RandLib::ur.irand(0,4)]; // pick a letter
//       int p =RandLib:: ur.irand(1,10); // pick power
//       p = (RandLib::ur.rand() > 0.7 ? p : -p);
//       // print
//       if (abs(p)>1)
// 	cout << l << "^{"<<p<<"} ";
//       else 
// 	cout << l << " ";
//     }
//     cout << "}$" <<endl << endl;
//   }

  
  char  alphabet[] = {'0','1'};

  
  for (int i=0;i<10;i++){ // make 10 sequences
    cout << "$\\mathbf{";
    for (int j=0;j<100;j++){ // of length 100
      char l = alphabet[RandLib::ur.irand(0,1)]; // pick a letter
      cout << l << " ";
    }
    cout << "}$" <<endl << endl;
  }
  return 0;
}
