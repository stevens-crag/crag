// Copyright (C) 2008 Alexander Ushakov
//
// Contents: Definition of an Image class
//
// Principal Author: Alexander Ushakov
//
// Status: in progress
//
// Revision History:
//


#ifndef _PDFgraphing_H_
#define _PDFgraphing_H_


#include "vector"
// #include "string"
#include "iostream"
using namespace std;

 
//---------------------------------------------------------------------------//
//--------------------------------- PDFPageObject ---------------------------//
//---------------------------------------------------------------------------//

 //! Implements interface for PDF drawing objects
class PDFPageObject
{
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
	
public:
  
  virtual void write( ostream& os ) { cout << "ERROR. Should not be used" << endl; }
};

//---------------------------------------------------------------------------//
//------------------------------- PDFPageObjectText -------------------------//
//---------------------------------------------------------------------------//


//! Class for pdf text object
class udPDFPageObjectText : public PDFPageObject
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:  

 //! Constructor of text pdf object
  /*! 
    \param x,y - coordinates where text will be placed
    \param t   - text to be output
    \param s   - text size (default 12)
  */
  udPDFPageObjectText( double x , double y , const string& t, int s = 12 ) : 
  	theX(x), theY(y), text(t), tSize(s) { }
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
	
public:
  
  void write( ostream& os );

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
    
 private:
  double theX, theY;
  string text;
  int tSize;
};


//---------------------------------------------------------------------------//
//------------------------------- PDFPageObjectLine -------------------------//
//---------------------------------------------------------------------------//

//! Class for pfd line object
class PDFPageObjectLine : public PDFPageObject
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:  
 //! Constructor  of a line pdf object
  /*! 
    \param x1,y1,x2,y2 - coordinates of the two endpoints of a line
  */  
  PDFPageObjectLine( double x1 , double y1 , double x2 , double y2 ) : theX1(x1), theX2(x2), theY1(y1), theY2(y2) { }
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
	
public:
  
  void write( ostream& os ) {
    os << theX1 << " " << theY1 << " m" << endl;
    os << theX2 << " " << theY2 << " l" << endl;
    os << "S" << endl;
  }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
    
 private:
  double theX1, theY1;
  double theX2, theY2;
};

//---------------------------------------------------------------------------//
//------------------------------- udPDFPageObjectLine -------------------------//
//---------------------------------------------------------------------------//

//! Class implements pdf line object
class udPDFPageObjectLine : public PDFPageObject
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

//! Constructor for pdf line object
  /*! 
    \param x1,y1,x2,y2 - coordinates of the two endpoints of a line
    \param gc          - gray scale intensity (0 - black, 1 - white)
    \param d           - if true than a dashed line is drawn
  */  
 public:  
  udPDFPageObjectLine( double x1 , double y1 , double x2 , double y2, double gc = 0, bool d = false ) : 
  theX1(x1), theX2(x2), theY1(y1), theY2(y2), gray_color( gc), dashed( d ) { }
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
	
public:
  
  void write( ostream& os );

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
    
 private:
  double theX1, theY1;
  double theX2, theY2;
  double gray_color;
  bool dashed;
};

//---------------------------------------------------------------------------//
//------------------------------- udPDFPageObjectVertLine -------------------//
//---------------------------------------------------------------------------//


//! Class implements a vertical line pdf object
class udPDFPageObjectVertLine : public udPDFPageObjectLine
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:  
 //! Constructor for vertical line pdf object
  /*! 
    \param x,y - coordinates of the top end point
    \param s   - length
    \param gc          - gray scale intensity (0 - black, 1 - white)
    \param d           - if true than a dashed line is drawn
  */  
 
  udPDFPageObjectVertLine( double x , double y , double s ,  double gc = 0, bool d = false ) :  
  			udPDFPageObjectLine( x , y , x , y + s, gc, d  )  { }
  
};

//---------------------------------------------------------------------------//
//------------------------------- udPDFPageObjectHorizLine -------------------//
//---------------------------------------------------------------------------//


//! Class implements a horizontal line pdf object
class udPDFPageObjectHorizLine : public udPDFPageObjectLine
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:  

 //! Constructor for horizontal line pdf object
  /*! 
    \param x,y - coordinates of the left end point
    \param s   - length
    \param gc  - gray scale intensity (0 - black, 1 - white)
    \param d   - if true than a dashed line is drawn
  */  
  udPDFPageObjectHorizLine( double x , double y , double s ,  double gc = 0, bool d = false ) :  
  			udPDFPageObjectLine( x , y , x + s , y, gc, d  )  { }
  
};

//---------------------------------------------------------------------------//
//------------------------------- udPDFPageObjectRect ---------------------//
//---------------------------------------------------------------------------//

//! Class implements a rectangle pdf object
class udPDFPageObjectRect : public PDFPageObject
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:  

 //! Constructor for a rectangle  pdf object
  /*! 
    \param x,y   - coordinates of the left end point
    \param w,h   - widht and height of the rectanlge
    \param bc  - border color from 0 to 1
    \param f   - fill color. If f < 0 then no filling
  */  
  udPDFPageObjectRect( double x , double y , double w , double h, double bc = 0, double f = -1 ) : 
  			theX(x), theY(y), theW(w), theH(h), border_color(bc), theFill( f )  { }
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
	
public:
  
  void write( ostream& os );

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
    
 private:
  double theX, theY;
  double theW, theH;
  double theFill;
  double border_color;
};

//---------------------------------------------------------------------------//
//------------------------------- udPDFPageObjectSquare ---------------------//
//---------------------------------------------------------------------------//


 //! Class implements a square pdf object
class udPDFPageObjectSquare : public udPDFPageObjectRect
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:  

 //! Constructor for a square pdf object
  /*! 
    \param x,y   - coordinates of the left end point
    \param s     - length of the sides
    \param bc    - border color from 0 to 1
    \param f     - fill color. If f < 0 then no filling
  */  
  udPDFPageObjectSquare( double x , double y , double s, double bc = 0,double f = -1 ) :
	  udPDFPageObjectRect( x , y , s , s, bc, f ) { }
  
};

//---------------------------------------------------------------------------//
//------------------------------- udPDFPageObjectCircle  --------------------//
//---------------------------------------------------------------------------//


 //! Class implements a square pdf object
class udPDFPageObjectCircle : public PDFPageObject
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:  

 //! Constructor for a circle pdf object
  /*! 
    \param cx,cy   - coordinates of the center of a circle
    \param r     - radius of a circle
    \param bc    - border color from 0 to 1
    \param f     - fill color. If f < 0 then no filling
  */  
  udPDFPageObjectCircle( double cx , double cy , double r, double bc = 0,double f = -1 ) : 
  theX(cx), theY(cy), theRad(r), border_color(bc), theFill( f )   { }
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
	
public:
  
  void write( ostream& os );

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
    
 private:
  double theX, theY;
  double theRad;
  double theFill;
  double border_color;
};

//---------------------------------------------------------------------------//
//------------------------------------ PDFPage ------------------------------//
//---------------------------------------------------------------------------//

//! Class implements a page of a pdf document
class PDFPage
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:  

  PDFPage( ) { }
  ~PDFPage( ) {
    for( vector< PDFPageObject* >::const_iterator o_it=theObjects.begin( ) ; o_it!=theObjects.end( ) ; ++o_it )
      delete *o_it;
  }

//! Write the contents of the page into a stream using PDF format 
  void writeContents( ostream& os ) {
    for( vector< PDFPageObject* >::const_iterator o_it=theObjects.begin( ) ; o_it!=theObjects.end( ) ; ++o_it )
      (*o_it)->write( os );
  }

//! Add a pdf object to a page
  void addObject( PDFPageObject* obj ) { theObjects.push_back( obj ); }
  
//! The actual  width of the page (no margins)
  static int width() { return 612; }
//! The actual  height of the page (no margins)
  static int height() { return 792; }

//! Width of the page excluding margins
  static int mWidth() { return width() - lMargin() - rMargin(); }
//! Height of the page excluding margins
  static int mHeight() { return height() - tMargin() - bMargin(); }
//! Coordinate of the right end of the page (excluding right margin)
  static int mRightEnd() { return mWidth() + lMargin(); }
//! Coordinate of the bottom end of the page (excluding the margin)
  static int mBottomEnd() { return mHeight() + tMargin(); }
  
//! Left margin  
  static int lMargin() { return 30; } 
//! Right margin  
  static int rMargin() { return 20; }
//! Top margin  
  static int tMargin() { return 20; } 
//! Bottom margin  
  static int bMargin() { return 20; }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
    
 private:
  vector< PDFPageObject* > theObjects;

};


//---------------------------------------------------------------------------//
//-------------------------------- PDFStructure -----------------------------//
//---------------------------------------------------------------------------//


//! Implements a pdf document consisting of several pages
class PDFStructure
{

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:  
  
  PDFStructure( ) { }
  ~PDFStructure( ) { 
    for( vector< PDFPage* >::iterator p_it=thePages.begin() ; p_it!=thePages.end( ) ; ++p_it )
      delete *p_it;
  }
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
	
public:

  //! Save pdf document into a file "filename"
  void save( const char* filename );
  //! Create a new page 
  void newPage( ) { 
    thePages.push_back( new PDFPage( ) );
  }

  //! Add a new object to the page
  /*! 
   \param p - page naumber
   \param obj - pointer to the new pdf object instance
   */
  void addObject( int p , PDFPageObject* obj ) { thePages[p]->addObject( obj ); }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions                                 //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:
  
  pair< const char* , int > preparePageContents( int p ) const;
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
    
 private:
 
  vector< PDFPage* > thePages;
  
};


#endif

