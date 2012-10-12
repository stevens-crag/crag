// Copyright (C) 2008 Alexey Myasnikov
// Contents: Example for class BraidDrawPDF
//
// Principal Authors: Alexey Myasnikov
//
// Revision History:
//

#include "WordDraw.h"
#include "Word.h"
#include <iostream>
#include "ShortBraidForm.h"
#include "RanlibCPP.h"

int main()
{
 
  int l = 1000;
  int N = 20;
  
  RandLib::ur.setSeedPID();
 
  // Construct a ramdom word in N-1 generators and compute its shortened braid form 
  Word w = Word::randomWord(N-1,l); 
  Word sw = shortenBraid(  N ,  w );
      
  //& Draw Braid Word ; How do I draw a braid word
  BraidDrawPDF bdpdf( N-1 );
  //& Draw Braid Word ; How do I change the size of cells in the table
  // Set the size of celss to 5 points
  bdpdf.setSS(5);

  // Compute the table for a word sw
  bdpdf.draw( sw );
  
  //& Draw Braid Word ; How do I save the image (in PDF format only)
  cout << "Saving to short_form.pdf ..." << endl;
  bdpdf.save("short_form.pdf");
   
  return 0;


}
