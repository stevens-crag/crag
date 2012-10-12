// Copyright (C) 2006 Alexei Miasnikov
//
// Contents: Definition of an Image class
//
// Principal Author: Alexei Miasnikov
//
// Status: in progress
//
// Revision History:
//


#ifndef _WORD_DRAW_H_
#define _WORD_DRAW_H_

#include "AImage.h"
#include "PDFgraphing.h"
#include "Word.h"
#include <fstream>
#include <string>

const int ssConst =  20;


//
//
//!  CREATES A PPM image of a table for braid word
//
//
class WordDraw
{
public:
  WordDraw( int n, const Word& w, bool draw_grid = true ): 
    ss( ssConst ),
    //    theWord(w),
    N(n),
    betBraids(0),
    theLength( 0 )
  {
    theImage = new CImage( w.length()*(ssConst+1) + 1 , n*(ssConst+1) + 1 );
    
    for (int i=0;i<theImage->getWidth();i++)
      for (int j=0;j<theImage->getHeight();j++){
	theImage->setRedPixel(i,j,255);
	theImage->setBluePixel(i,j,255);
	theImage->setGreenPixel(i,j,255);
      }
    
    drawCompressedBraid( w );

    if (draw_grid){
      for (int i=0;i<=(theLength+1)*(ss+1);i+=(ss+1))
	drawVerticalGrid( i, (i/(ss+1) % 10 ? 220 :0) );
      for (int i=0;i<theImage->getHeight();i+=(ss+1))
	drawHorizontalGrid( i );
    }
  
  }

  WordDraw( int n, const list<Word>& w, bool draw_grid = true ): 
    ss( ssConst ),
    //    theWords(w),
    N(n),
    betBraids(20),
    theLength( 0 )
  {
    int max_length = 0;
    for (list<Word>::const_iterator I=w.begin();I!=w.end();I++)
      if (max_length  < I->length())
	max_length = I->length();
    
    int one_braid_width = n*(6+1)+1;
    theImage = new CImage( max_length*(ssConst+1) + 1 , w.size()*one_braid_width + 1 + (w.size()-1)*betBraids );
      
    for (int i=0;i<theImage->getWidth();i++)
      for (int j=0;j<theImage->getHeight();j++){
	theImage->setRedPixel(i,j,255);
	theImage->setBluePixel(i,j,255);
	theImage->setGreenPixel(i,j,255);
      }
    
    int i_offset=0;
    for (list<Word>::const_iterator I=w.begin();I!=w.end();I++,i_offset+=one_braid_width+betBraids)
      drawCompressedBraid( *I, i_offset );

    if (draw_grid){
      for (int i=0;i<=(theLength+1)*(ss+1);i+=(ss+1))
	drawVerticalGrid( i , ( i/(ss+1) % 10 ? 220 : 0 ));
      
      i_offset=one_braid_width;
      for (int i=0;i<theImage->getHeight();i+=(ss+1)){
	if ( i >= i_offset){
	  i+=betBraids-ss;
	  i_offset+=one_braid_width+betBraids;
	}
	drawHorizontalGrid( i );
      }
    }
  
  }
  ~WordDraw() { delete theImage; }
  void saveTo( const string& f_name ) {theImage->saveTo(f_name);}
 
private:
  //  void drawBraid(){
    
  //    int i=0;
  //    for (ConstWordIterator wI=theWord.begin();wI!=theWord.end();wI++,i++)
  //      drawGenerator( *wI,i );    
  //  }
  
  
  void drawCompressedBraid( const Word& theWord, int vert_offset = 0){
    vector<int> positions(N+2,-1); // +2 for the left and right padding
    
    int i=0;
    for (ConstWordIterator wI=theWord.begin();wI!=theWord.end();wI++,i++){
      Generator g = *wI;
      int absg = abs(g);
      int pos = max(positions[absg],max(positions[absg-1],positions[absg+1])) + 1;
      drawGenerator( g,pos,vert_offset );
      positions[absg] = pos;
      if (theLength < pos) theLength = pos;
    }
    
  }
  
  void drawGenerator( Generator g, int pos, int vert_offset ){
    for (int iss=0;iss<ss;iss++)
      for (int jss=0;jss<ss;jss++)
	if ( g > 0 ){
	  theImage->setRedPixel( pos*ss + iss + pos+1,  vert_offset+ g + (g - 1)*ss  + jss,168);
	  theImage->setBluePixel( pos*ss + iss + pos+1, vert_offset+ g + (g - 1)*ss  + jss,168);
	  theImage->setGreenPixel( pos*ss + iss + pos+1, vert_offset+g + (g - 1)*ss  + jss,168);
	} else {
	  int ag = abs(g);
	  theImage->setRedPixel( pos*ss +iss + pos+1,  vert_offset+ ag + (ag-1)*ss + jss,0); 
	  theImage->setBluePixel( pos*ss +iss + pos+1, vert_offset+ ag + (ag-1)*ss + jss,0); 
	  theImage->setGreenPixel( pos*ss +iss + pos+1,vert_offset+ ag + (ag-1)*ss + jss,0); 
	}
  }
  
  void drawVerticalGrid( int vpos, int color){ 
    
    for (int i=0;i<theImage->getHeight();i+=2){
   	  theImage->setRedPixel( vpos, i,color); 
	  theImage->setBluePixel( vpos, i, color );
	  theImage->setGreenPixel( vpos, i, color );
    }
  }
  
  void drawHorizontalGrid( int hpos){
    for (int i=0;i<=(theLength+1)*(ss+1);i+=2){
   	  theImage->setRedPixel( i, hpos,220); 
	  theImage->setBluePixel( i, hpos, 220 );
	  theImage->setGreenPixel( i, hpos, 220 );
    }
  }
  
  //  const Word& theWord;
  int N;
  CImage* theImage;
  int ss;
  int theLength;
  int betBraids;
};

//
//
//  CREATES A PS image of a table for braid word
//
//

class RGB
{
public:
	RGB(double r, double g, double  b): R(r), G(g), B(b) {} 
	RGB(double gray ): R(gray),G(gray),B(gray) {}


	double R;
	double G;
	double B; 
};


//! Class for drawing a braid on n strands as a square table 
class BraidDrawPDF
{
public:

//! Constructor
/*!
 \param n - number of generators (i.e. number of strands -1 )
 \param t - optional title which will appear on each page 
 */
  BraidDrawPDF( int n, string t = string() ):
    ss( ssConst ),
    //    theWord(w),
    N(n),
    betBraids(0),
    theLength( 0 ),
    theTitle( t ),
    useCircle( false ) // @am There is a problem with using circles when words size is large (say 10000)
  	{
  	}
  	
  	~BraidDrawPDF() {  }

//! Draw a braid represented by the word w
void draw( const Word& w) {
	// need to start a new page!!!!
	pdf_out.newPage();
   	drawCompressedBraid( w );

}


//! Write the picture into a file
void save( const char* f ){
 	pdf_out.save(f);
 }
	
//! Set the size of the square cells of the table
 void setSS(int new_ss ) { ss = new_ss; }


private:

 void drawGrid( int nP, int sN) {
   	for (int i=PDFPage::lMargin();i<=PDFPage::lMargin()+PDFPage::mWidth();i+=(ss)) // )i<=(posInPage()+1)*(ss);i+=(ss)) // theLength
		pdf_out.addObject( nP , new udPDFPageObjectVertLine( i, stripPosition( sN ),tableWidth(), 0.9, true ) );
    for (int i=stripPosition( sN );i<stripPosition( sN + 1);i+=(ss))
		pdf_out.addObject( nP , new udPDFPageObjectHorizLine( PDFPage::lMargin() ,i,PDFPage::mWidth(), 0.9, true ) );		
	
	// Draw a line after the grid
	pdf_out.addObject( nP , new udPDFPageObjectRect(PDFPage::lMargin() , stripPosition( sN + 1) - ss, PDFPage::mWidth(), ss, 0.95, 0.95 ) );
    
    // add grid mark
//    pdf_out.addObject( nP , new PDFPageObjectText( 10 , stripPosition( sN ) , string("AAAAA") ) );
    pdf_out.addObject( nP , new udPDFPageObjectText( PDFPage::lMargin() , PDFPage::tMargin() + 10 , theTitle ) );

  }
  
 int tableWidth() const { return N*(ss);}
 //int stripPosition(int i){return PDFPage::height() - i*tableWidth() - i*2;}
 int stripPosition(int i){return i*tableWidth() + i*ss + PDFPage::tMargin() + (theTitle.length() > 0 ? 20 : 0) ;} // i*ss is the gap between stripes + top margin
 
 
 
 
private:
  
  void drawCompressedBraid( const Word& theWord, int vert_offset = 0){
    vector<int> positions(N+2,-1); // +2 for the left and right padding
    
 //   drawGrid(0,0);
	
	cout << "Positions in page: " << posInPage() << endl;
    cout << "Stripes per page : " << stripesInPage() << endl;
    
    int i=0;
    int maxSN = 0;
    int maxPage = 0;
    for (ConstWordIterator wI=theWord.begin();wI!=theWord.end();wI++,i++){
   	  Generator g = *wI;
      int absg = abs(g);
      int pos = max(positions[absg],max(positions[absg-1],positions[absg+1])) + 1;
      
          
      if ( maxSN < int(pos / posInPage()) ){   	
  		drawGrid(maxPage,maxSN % stripesInPage() );
  	  	maxSN = int(pos / posInPage());
//     	drawGrid(maxPage,maxSN % stripesInPage() );
      }
      if ( maxPage < maxSN / stripesInPage() ){
  		pdf_out.newPage();    	
      	maxPage++;
//     	drawGrid(maxPage,maxSN % stripesInPage() );
     	//     	cout << "MAX PAGE " << maxPage << endl;
      }
      
      positions[absg] = pos;
      if (theLength < pos) theLength = pos;    
    
     drawGenerator( g,pos,vert_offset );
   
 
    }
    drawGrid(maxPage,maxSN % stripesInPage() );   
  }

  int posInPage() const { return PDFPage::mWidth() / (ss); }
  int stripesInPage() const { return PDFPage::mHeight() / ((N+1)*(ss)); }
  	
  void drawGenerator( Generator g, int pos, int vert_offset ){
  	
  	int sPos  = pos % posInPage(); // position in the strip 
 
	int pN = getPageStrip( pos ).first;
	int sN = getPageStrip( pos ).second;
 	int new_vert_offset = stripPosition(sN);
  	
	if ( g > 0 ){
		if (useCircle)
			pdf_out.addObject( pN , new udPDFPageObjectCircle(PDFPage::lMargin() + sPos*ss + double(ss)/2.0 , new_vert_offset + (g-1)*(ss) + double(ss)/2.0, double(ss)/2.0, 0.9, 0.5 ) );
		else
			pdf_out.addObject( pN , new udPDFPageObjectSquare(PDFPage::lMargin() + sPos*ss , new_vert_offset + (g-1)*(ss), ss, 0.9, 0.5 ) );

	} else {
	  	int ag = abs(g);
		if (useCircle)
			pdf_out.addObject( pN , new udPDFPageObjectCircle(PDFPage::lMargin() + sPos*ss + double(ss)/2.0, new_vert_offset + (ag-1)*(ss)+double(ss)/2.0, double(ss)/2.0, 0.9, 0 ) );
		else
			pdf_out.addObject( pN , new udPDFPageObjectSquare(PDFPage::lMargin() + sPos*ss, new_vert_offset + (ag-1)*(ss), ss, 0.9, 0 ) );
	}
  }
  
  pair<int,int> getPageStrip( int pos ) const {
 	int totSN = pos / posInPage();  // the strip number 	
  	int pN = totSN / stripesInPage(); // the page number
	int sN = totSN % stripesInPage(); // the strip number on page pN
 	
//	cout << pN << endl;
 	
 	return pair<int,int>(pN,sN);
  } 

  //  const Word& theWord;
  int N;
  PDFStructure pdf_out;
  int ss;
  int theLength;
  int betBraids;
  bool useCircle;
  string theTitle;
};
#endif
