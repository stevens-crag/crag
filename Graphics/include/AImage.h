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


#ifndef _IMAGE_H_
#define _IMAGE_H_

#include <iostream>
#include <string>

using namespace std;
const int MAXGRAY=256;


enum IMAGE_TYPE { GRAY, COLOR };
enum FILE_TYPE { BW, PGM, PPM, JPG, BMP, FT_NA };


FILE_TYPE getFileType(const string& in_file_name);

// --------------- LookUpTable class --------------------- //

class LookUpTable
{
 public:
  LookUpTable();
  // default constructor

  LookUpTable( const LookUpTable&);
  // Copy constructor

  LookUpTable& operator = (const LookUpTable&);

  void set(int,int);
  
  int get(int) const;
 private:
  unsigned char table[256];
};

// --------------- Abstract Image class ----------------------------- //

class AImage
{
 public:

  AImage():  width(0),  height(0), size(0){ }
  // Default constructor creates empty image
  
  AImage(int w,int h)
	  :   width(w),  height(h), size( w*h) { }
  // Create an empty image of w x h size

  AImage( const  string& in_file_name) : width(0),  height(0), size(0){ }
  // Create image from the input file
  
  AImage( const AImage& i) : width(i.width),  height(i.height), size(i.width * i.height) { }
  // Copy constructor
  
  AImage( const AImage& image, int tlx, int tly, int brx, int bry, bool zero_pading) 
	  : width(image.width + tlx + brx),
		height(image.height + bry + bry),
		size((image.width + tlx + brx) * (image.height + tly + bry))
{ }
  // Copy image with padding

  int getSize() const { return size; }

  int getWidth() const { return width; }

  int getHeight() const { return height; }

 // virtual void printOn( CArchive& out) const;

  virtual void saveTo( const string& file_name ) = 0;
  
  virtual IMAGE_TYPE getType() const = 0;

protected:

  int width;
  // The width of the image
  int height;
  // The height of the image
  int size;
  // image size
  
};


// ----------------------- LUTImage class ------------------- //

class LUTImage 
{
public:

  LUTImage() { }
  // Default constructor creates empty image
  
  LUTImage(const AImage& i) { }
  // Default constructor creates empty image
  
  LUTImage(int w, int h) { }
  
  LUTImage(const LUTImage& li) : theLookUpTable( li.theLookUpTable ) { }

  void setLookUpTable( const LookUpTable& t){
    theLookUpTable = t;
  }

  void setInLT( int i, int v){
    theLookUpTable.set(i,v);
  }

  int getInLT( int i) const {
    return theLookUpTable.get(i);
  }

protected:

  LookUpTable theLookUpTable;
};


// --------------- Image class ----------------------------- //

class GRImage : public AImage
{
 public:

  GRImage(): AImage(), image_string( NULL ) { }
  // Default constructor creates empty image
  
  GRImage(int w,int h): AImage(w,h), image_string( new unsigned char[h*w] ) 
    {
      for (int i=0;i<size;i++)
	image_string[i] = 0;
    }
  // Create an empty image of w x h size

  GRImage( const string& in_file_name);
  // Create image from the input file
  
  GRImage( const GRImage& image);
  // Copy constructor
  
  GRImage( const GRImage&, int, int , int , int , bool );
  // Copy image with padding

  ~GRImage( );

  void setPixel(int i, int j,int v);

  virtual int getPixel(int i, int j) const;

  void setPixel(int n,int v);

  void setPixelCliped(int n,int v);

  void setPixelCliped(int i,int j,int v);

  virtual int getPixel(int n) const;

 // virtual void printOn( CArchive& out) const;

  void saveTo( const string& file_name );

  IMAGE_TYPE getType() const {return GRAY; }
  
  GRImage& operator = (const GRImage&);

protected:

  friend class GRLUTImage;	 
  
  virtual void printOnPGM( ostream& out)const;
  virtual void printOnBW( ostream& )const;
  // Output the image into a stream
  
  unsigned char* image_string;
  // Image data 

  void readFromBW( istream& in);
  void readFromPGM( istream& in);
};


// ----------------------- LUTImage class ------------------- //

class GRLUTImage : public GRImage, public LUTImage
{
public:
  
  
  GRLUTImage(const GRLUTImage& li) : GRImage( li ), LUTImage((const LUTImage&) li ) { 	}
  
  
  GRLUTImage() : GRImage() { }
  // Default constructor creates empty image
  
  GRLUTImage(const GRImage& i) : GRImage( i ) { }
  // Default constructor creates empty image
  
  GRLUTImage(int w, int h) : GRImage(w,h) { }
  // Create an empty image of w x h size
  
  int getPixel(int i, int j) const;
  
  int getPixel(int n) const;
  
  
  inline friend ostream& operator << ( ostream& o, const GRLUTImage& i ) {
    i.printOn( o );
    return o;
  }
  // Output operator - prints data in ostream


private:

  void printOn( ostream& )const;
};




// --------------- Color Image class ----------------------------- //

class CImage : public AImage
{
 public:
  CImage(): AImage() { }
  // Default constructor creates empty image
  
  CImage(int w,int h): AImage(w,h), redImage(w,h), blueImage(w,h), greenImage(w,h) { }
  // Create an empty image of w x h size

  CImage( const string& in_file_name, FILE_TYPE ft = PPM);
  // Create image from the input file
  
  CImage( const CImage& image) : 
    AImage(image.width, image.height),
    redImage( image.redImage ),
    blueImage( image.blueImage ),
    greenImage( image.greenImage ) { }
  // Copy constructor
  
  CImage( const CImage& i, int p1, int p2, int p3, int p4, bool p5) : 
    AImage(i, p1, p2, p3, p4, p5),
    redImage(i.redImage, p1, p2, p3, p4, p5),
    blueImage(i.blueImage, p1, p2, p3, p4, p5),
    greenImage(i.greenImage, p1, p2, p3, p4, p5) { }
  // Copy image with padding
  
  // modifiers
  /////////////////////////////////////
  
  void setRedPixel(int i, int j,int v) { redImage.setPixel(i,j,v); }
  void setRedPixel(int n,int v) { redImage.setPixel(n,v); }

  void setBluePixel(int i, int j,int v) { blueImage.setPixel(i,j,v); }
  void setBluePixel(int n,int v) { blueImage.setPixel(n,v); }
  
  void setGreenPixel(int i, int j,int v) { greenImage.setPixel(i,j,v); }
  void setGreenPixel(int n,int v) { greenImage.setPixel(n,v); }

  void setRedPixelCliped(int n,int v) { redImage.setPixelCliped(n,v); }
  void setRedPixelCliped(int i,int j,int v) { redImage.setPixelCliped(i,j,v); }

  void setBluePixelCliped(int n,int v) { blueImage.setPixelCliped(n,v); }
  void setBluePixelCliped(int i,int j,int v){ blueImage.setPixelCliped(i,j,v); }

  void setGreenPixelCliped(int n,int v) { greenImage.setPixelCliped(n,v); }
  void setGreenPixelCliped(int i,int j,int v) { greenImage.setPixelCliped(i,j,v); }


  // accessors
  virtual unsigned char getRedPixel(int i, int j) const { return redImage.getPixel(i,j); }
  virtual unsigned char getRedPixel(int n) const { return redImage.getPixel(n); }

  virtual unsigned char getBluePixel(int i, int j) const { return blueImage.getPixel(i,j); }
  virtual unsigned char getBluePixel(int n) const { return blueImage.getPixel(n); }

  virtual unsigned char getGreenPixel(int i, int j) const { return  greenImage.getPixel(i,j); }
  virtual unsigned char getGreenPixel(int n) const {  return greenImage.getPixel(n); }

  const GRImage& getRedImage() const { return redImage; }

  const GRImage& getGreenImage() const { return greenImage; }

  const GRImage& getBlueImage() const { return blueImage; }

  // virtual void printOn( CArchive& out) const;

  void saveTo( const string& file_name );

  IMAGE_TYPE getType() const { return COLOR; }
  

protected:

  friend class CLUTImage;	 
  CImage& operator = (const CImage&);
  // Assign operator is prohibited
  
  // Output the image into a stream
  void printOnPPM( ostream& out)const;

  void readFromPPM( istream& in);

  // Image data 
  GRImage redImage;
  GRImage blueImage;
  GRImage greenImage;
};


// ----------------------- Color LUTImage class ------------------- //

class CLUTImage : public CImage, public LUTImage
{
public:

	CLUTImage(const CLUTImage& li) : CImage( li ), LUTImage((const LUTImage&)li){ } 
	
	CLUTImage() : CImage() { }
  // Default constructor creates empty image
  
	CLUTImage(const CImage& i) : CImage( i ){ }
	

	CLUTImage(int w, int h) : CImage(w,h) { }
  // Create an empty image of w x h size


	//modifiers
 
  // modifiers
  /////////////////////////////////////

 
  // accessors
  virtual unsigned char getRedPixel(int i, int j) const { 
	  return   theLookUpTable.get(CImage::getRedPixel(i,j));
  }
  virtual unsigned char getRedPixel(int n) const { 
	  return   theLookUpTable.get(CImage::getRedPixel(n));
  }

  virtual unsigned char getBluePixel(int i, int j) const { 
	  return   theLookUpTable.get(CImage::getBluePixel(i,j));
  }
  virtual unsigned char getBluePixel(int n) const { 
	  return   theLookUpTable.get(CImage::getBluePixel(n));
  }
  virtual unsigned char getGreenPixel(int i, int j) const { 
	  return   theLookUpTable.get(CImage::getGreenPixel(i,j));
  }
  virtual unsigned char getGreenPixel(int n) const{ 
	  return   theLookUpTable.get(CImage::getGreenPixel(n));
  } 


private:
  
//  void printOnPPM(ostream& out)const;
};





void convert(const CImage* ci, GRImage* gi);


#endif
