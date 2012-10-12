// Copyright (C) 2008 Alexander Ushakov
// Contents: Implementation of class 
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "PDFgraphing.h"
#include "strstream" 
#include "iostream" 
#include "fstream" 
#include "vector" 
#include "map" 
#include <math.h>

using namespace std;


const double LINE_WIDTH = 0;

//---------------------------------------------------------------------------//
//------------------------------- udPDFPageObjectText -----------------------//
//---------------------------------------------------------------------------//

void udPDFPageObjectText::write( ostream& os ) {
	os << 0 << " g" << endl;
  	os << "BT" << endl
	   << "/F1 " << tSize << " Tf" << endl
	   << theX << " " << PDFPage::height()  - theY << " Td" << endl
	   << "( " << text << " ) Tj" << endl
	   << "ET" << endl;
  }
  
  
//---------------------------------------------------------------------------//
//------------------------------- udPDFPageObjectLine -----------------------//
//---------------------------------------------------------------------------//


void udPDFPageObjectLine::write( ostream& os )
 {
 	os << LINE_WIDTH << " w" << endl;
 	if (dashed)
//	 	os << "[2 3] 11 d" << endl; 
	 	os << "[2] 1 d" << endl; 
 	os << gray_color << " G" << endl;
 	os << theX1 << " " << PDFPage::height() - theY1 << " m" << endl;
    os << theX2 << " " << PDFPage::height() - theY2 << " l" << endl;
    os << "S" << endl;
 	if (dashed)
	 	os << "[] 0 d" << endl; 
    
  }

//---------------------------------------------------------------------------//
//------------------------------- udPDFPageObjectRect -------------------------//
//---------------------------------------------------------------------------//

void udPDFPageObjectRect::write( ostream& os ) {

	os  << theX << " " << PDFPage::height() - (theY + theH) << " "
		<< theW << " " <<  theH << " "
		<< "re" << endl;
	
	os << LINE_WIDTH << " w" << endl;
 	os << border_color << " G" << endl;

	if (theFill < 0 ) {
	   os << "s" << endl;
	} else {
	   
	   os << theFill << " g" << endl
	      << "b" << endl;
	      
	}
  }
  
//---------------------------------------------------------------------------//
//-------------------------------- PDFStructure -----------------------------//
//---------------------------------------------------------------------------//


pair< const char* , int > PDFStructure::preparePageContents( int p ) const
{
  const int buffer_size = 100*4096;
  static char buffer[buffer_size];

  ostrstream ostr( buffer , buffer_size , ios::out );
  thePages[p]->writeContents( ostr );
  return pair< const char* , int >( buffer , ostr.tellp( ) );
}

//---------------------------------------------------------------------------//
//------------------------------- udPDFPageObjectCirce -------------------------//
//---------------------------------------------------------------------------//

void udPDFPageObjectCircle::write( ostream& os ) {

	os << LINE_WIDTH << " w" << endl;
 	os << border_color << " G" << endl;

    os << theX+theRad << " " << PDFPage::height() - theY << " m" << endl;
    
    // diamond
//    os << theX << " " << theY - theRad<< " l" << endl;
//    os << theX - theRad << " " << theY  << " l" << endl;
//    os << theX  << " " << theY + theRad << " l" << endl;
    
    // do circle path
    for (double angle = 0.1; angle < 2*M_PI;angle+=0.5){
    
    	os << theX + cos(angle)*theRad<< " " << PDFPage::height() - (theY + sin(angle)*theRad)<< " l" << endl;
    	
    }
    
	if (theFill < 0 ) {
	   os << "s" << endl;
	} else {
	   
	   os << theFill << " g" << endl
	      << "b" << endl;
	      
	}
  }



//---------------------------------------------------------------------------//
//-------------------------------- PDFStructure -----------------------------//
//---------------------------------------------------------------------------//


void PDFStructure::save( const char* filename )
{
  // Aux buffer, string-stream, and the file
  int buffer_size = 100*4096;
  char buffer[buffer_size];
  // char buffer2[buffer_size];
  // ostrstream ostr2( buffer2 , buffer_size );
  ofstream OF( filename );

  
  // Object offsets
  int curOffset = 0;
  vector< int > theOffsets;

  
  // Header
  OF << "\%PDF-1.4" << endl;
  theOffsets.push_back( OF.tellp( ) );
  
  // Write the catalog entry
  OF << "1 0 obj" << endl;
  OF << "<< /Type /Catalog" << endl;
  OF << "/Outlines 2 0 R" << endl;
  OF << "/Pages 3 0 R" << endl;
  OF << ">>" << endl;
  OF << "endobj" << endl;
  theOffsets.push_back( OF.tellp( ) );
  
  // Write the Outline entry
  OF << "2 0 obj" << endl;
  OF << "<< /Type /Outlines" << endl;
  OF << "/Count 0" << endl;
  OF << ">>" << endl;
  OF << "endobj" << endl;
  theOffsets.push_back( OF.tellp( ) );

  // the number of objects in the header, i.e. before pages 
  int headObjectNum = 5; 


  // Write the description of pages
  OF << "3 0 obj" << endl;
  OF << "<< /Type /Pages" << endl;
  OF << "/Kids [" << endl;
  for( int p=0 ; p<thePages.size( ) ; ++p )
    OF << headObjectNum+p << " 0 R" << endl;
  OF << "]" << endl;
  OF << "/Count " << thePages.size( ) << endl;
  OF << ">>" << endl;
  OF << "endobj" << endl << endl;
  theOffsets.push_back( OF.tellp( ) );
  
  // Write Font entry
  OF << "4 0 obj" << endl
	 << "<< /Type /Font" << endl
	 << "/Subtype /Type1" << endl
	 << "/Name /F1" << endl
	 << "/BaseFont /Helvetica" << endl
	 << "/Encoding /MacRomanEncoding" << endl
	 << ">>" << endl
	 << "endobj" << endl << endl;
  theOffsets.push_back( OF.tellp( ) );
  
  
  // Write an object for each page
  for( int p=0 ; p<thePages.size( ) ; ++p ) {
    // OF << "% ====== PAGE N " << p << " =====" << endl;
    OF << headObjectNum+p << " 0 obj" << endl;
    OF << "<< /Type /Page" << endl;
    OF << "/Parent 3 0 R" << endl;
    OF << "/MediaBox [0 0 612 792]" << endl;
    OF << "/Contents " << headObjectNum+thePages.size( )+p << " 0 R" << endl;
    OF << "/Resources << /ProcSet [/PDF /Text] " << endl;
    OF << "/Font << /F1 4 0 R >>" << endl << ">>" << endl;
    OF << ">>" << endl;
    OF << "endobj" << endl << endl;
    theOffsets.push_back( OF.tellp( ) );
  }
  
  // Write the contents for each page
  for( int p=0 ; p<thePages.size( ) ; ++p ) {
    // OF << "% ====== PAGE N " << p << " contents =====" << endl;
    OF << headObjectNum+thePages.size( )+p << " 0 obj" << endl;
    pair< const char* , int > contents = PDFStructure::preparePageContents( p );
    OF << "<< /Length " << contents.second << " >>" << endl;
    OF << "stream" << endl;
    OF.write( contents.first , contents.second );
    OF << "endstream" << endl;
    OF << "endobj" << endl << endl;
    theOffsets.push_back( OF.tellp( ) );
  }
  

  // Write Cross-reference table
  OF << "xref" << endl;
  OF << "0 " << theOffsets.size() << endl;
  OF << "0000000000 65535 f" << endl;
  OF.fill( '0' );
  for( int i=0 ; i<theOffsets.size()-1 ; ++i ) {
    OF.width( 10 );
    OF << theOffsets[i];
    OF << " 00000 n" << endl;
  }
  OF << "trailer" << endl;
  OF << "<< /Size " << theOffsets.size() << endl;
  OF << "/Root 1 0 R" << endl;
  OF << ">>" << endl;
  OF << "startxref" << endl;
  OF << theOffsets[theOffsets.size()-1] << endl;
  
  
  // end of file
  OF << "\%\%EOF" << endl;
}



