// Copyright (C) 2006 Alexei Miasnikov
//
// Contents: Implementation of Image classes
//
// Principal Author: Alexei Miasnikov
//
// Status: in progress
//
// Revision History:
//


#include "AImage.h"
#include <string.h>
#include "errormsgs.h"
#include <fstream>
#include <sstream>

FILE_TYPE getFileType(const string& in_file_name)
{
  int pl = in_file_name.length() - 1;
  char ext[4];
  int	count = pl;
  while (in_file_name[count] != '.' && count > 0){
    ext[pl - count] = in_file_name[count];
    count--;
  }
  ext[pl - count] = '\0';
  
  if (strcmp(ext,"mpp") == 0 || strcmp(ext,"MPP") == 0 )
    return PPM;
  else if ( strcmp(ext,"wb") == 0 || strcmp(ext,"WB") == 0)
    return BW;
  else if ( strcmp(ext,"mgp") == 0 || strcmp(ext,"MGP")==0  )
    return PGM;
  else if ( strcmp(ext,"pmb") == 0 ||  strcmp(ext,"PMB") == 0  )
    return BMP;
  else if ( strcmp(ext,"gpj") == 0 || strcmp(ext,"GPJ") ==0  )
    return JPG;
  else
    return FT_NA;
}


//  Class implemetnation

LookUpTable::LookUpTable()
  // default constructor
{
  for (int i=0;i<256;i++)
    table[i] = (unsigned char)i;
}

LookUpTable::LookUpTable( const LookUpTable& t)
  // Copy constructor
{
  for (int i=0;i<256;i++)
    table[i] = t.table[i];
}

LookUpTable& LookUpTable::operator = (const LookUpTable& t)
{
  for (int i=0;i<256;i++)
    table[i] = t.table[i];
  return *this;
}


void LookUpTable::set(int i,int v)
{
  if ( i < 0 || i >= 256){ 
    msgs::error("void LookUpTable::set(char i,char v): index out of bounds");
  }
  
  if ( v <= 0 )
    table[i] = 0;
  else if ( v >= 255 )
    table[i] = 255;
  else
    table[i] = v;
  
}

int LookUpTable::get(int i) const
{
  // error_check( i < 0 || i >= 256, 
  //	       "void LookUpTable::set(char i,char v): "
  //	       "index out of bounds");
  return table[i];
}

//////////////////////////////////////////////////////////////////////
//
//  GRIMAGE
//
/////////////////////////////////////////////////////////////////////


GRImage::GRImage( const string& in_file_name):AImage( in_file_name )
  // Create image from the input file
{
  ifstream in(in_file_name.c_str(),ios::in |ios::binary);
  int pl = in_file_name.length() - 1;
  char ext[4];
  int count = pl;
  while (in_file_name[count] != '.'){
    ext[pl - count] = in_file_name[count];
    count--;
  }
  ext[pl - count] = '\0';
  
  if (strcmp(ext,"wb") == 0 )
    readFromBW(in);
  else if(strcmp(ext,"mgp") == 0 )
    readFromPGM(in);
  in.close();
}


GRImage::GRImage( const GRImage& image) :AImage( image )
  // Copy constructor
{
  image_string = new unsigned char[size];
  for (int i=0;i<size;i++){
    image_string[i] = image.getPixel(i);
  }
}

GRImage::~GRImage( )
{
  if (image_string)
    delete [] image_string;
}


GRImage::GRImage( const GRImage& image, int tlx, int tly, int brx, int bry, bool zero_pading)
	:AImage(image, tlx, tly, brx, bry, zero_pading)
{
	image_string = new unsigned char[size];
	// set image in the center
	for (int j=tlx;j<width-brx;j++){
		for (int i=tly;i<height-bry;i++)
			setPixel(j,i,image.getPixel(j-tlx,i-tly));
	}
	
	// install padding
	// collomns on right and left
	for (int h1=tly;h1<height-bry;h1++){
		for (int i1=0;i1<tlx;i1++)
			setPixel(i1,h1,(zero_pading ? 0 : getPixel(tlx,h1)) );
		for (int i2=width-brx-1;i2<width;i2++)
			setPixel(i2,h1,(zero_pading ? 0 : getPixel(width-brx-1,h1)) );
	}
	// top and bottom rows
	for (int w1=0;w1<width;w1++){
		for (int i1=0;i1<tly;i1++)
			setPixel(w1,i1,(zero_pading ? 0 : getPixel(w1,tly)) );
		for (int i2=height-bry-1;i2<height;i2++)
			setPixel(w1,i2,(zero_pading ? 0 : getPixel(w1,height-bry-1)) );
	}
}


GRImage& GRImage::operator = (const GRImage& i)
{
	if (image_string != NULL )
		delete [] image_string;

	width = i.width;
	height = i.height;
	size = i.size;

	image_string = new unsigned char[size];
	for (int ii=0;ii<size;ii++){
	    image_string[ii] = i.getPixel(ii);
	}
	return *this;
}


void GRImage::setPixel(int i, int j, int v){
  /*  error_check(i < 0 || i>=width,"Image::setPixel(int i, int j) const : "
	      "index out of bounds.");
  error_check(j < 0 || j>=height,"Image::setPixel(int i, int j) const : "
	      "index out of bounds.");
  error_check( v < 0 || v >=256,"Image::setPixel(int i, int j) const : "
	       "intensity out of bounds.");
  */
  int n = j*width+i;
  image_string[n] = v;
}

void GRImage::setPixel(int n, int v)
{
  /*
  error_check(n >=size || n < 0,"Image::setPixel(int i, int j) const : "
	      "index out of bounds.");
  error_check( v < 0 || v >=256,"Image::setPixel(int i, int j) const : "
	       "intensity out of bounds.");
 */
  image_string[n] = v;
}

void GRImage::setPixelCliped(int n, int v)
{
  int intensity = v;
  if ( intensity < 0 )
    intensity = 0;

  if ( intensity > 255 )
    intensity = 255;
  
  setPixel( n, intensity);
}

void GRImage::setPixelCliped(int i,int j, int v)
{
  int n = j*width+i;
  setPixelCliped(n, v);
}

int GRImage::getPixel(int i, int j) const
{
  /*  error_check(i < 0 || i>=width,"Image::setPixel(int i, int j) const : "
	      "index out of bounds.");
  error_check(j < 0 || j>=height,"Image::setPixel(int i, int j) const : "
	      "index out of bounds.");
  */
  int n = j*width+i;
  
  //  error_check(n >=size || n < 0,"Image::setPixel(int i, int j) const : "
  //      "index out of bounds."); 
  return image_string[n]; 
}

int GRImage::getPixel(int n) const
{
  //  error_check(n >=size || n < 0,"Image::setPixel(int i, int j) const : "
  //	      "index out of bounds."); 
  return image_string[n]; 
}

void GRImage::printOnPGM( ostream& out)const
  // Output the image into a stream
{
  out << "P5" << endl << "#MAD Alex" << endl<< width <<" "<<height<<endl <<255 << endl;
  
  out.flush();
  out.write((char*)image_string,size);
  out.flush(); 
}


void GRImage::printOnBW( ostream& out)const
  // Output the image into a stream
{
  
  int w = width;
  int h = height;
  unsigned char p[8];
   
  p[0] = ((w>>24) & 0xff);	/* save 4-byte width  */
  p[1] = ((w>>16) & 0xff);
  p[2] = ((w>> 8) & 0xff);
  p[3] = ( w      & 0xff);
  
  p[4] = ((h>>24) & 0xff);	/* save 4-byte height */
  p[5] = ((h>>16) & 0xff);
  p[6] = ((h>> 8) & 0xff);
  p[7] = ( h      & 0xff);

  out.write((char*)p,8);
  out.flush();
  out.write((char*)image_string,size);
  out.flush(); 
}


void GRImage::readFromBW(istream &in)
{
  unsigned char p[8];
  in.read((char*)p, 8);

  width = (p[0]<<24) + (p[1]<<16) + (p[2]<<8) + p[3];
  height = (p[4]<<24) + (p[5]<<16) + (p[6]<<8) + p[7];
  
  size = width*height;
  image_string = new unsigned char[size];

  in.read((char*)image_string,size);

}

void GRImage::readFromPGM(istream &in)
{
	char tmpLine[1000];
	in.getline(tmpLine,1000);
	if (strcmp(tmpLine, "P5") != 0)
	  msgs::error("Can't understand the header");
	tmpLine[0] = '#';
	while( tmpLine[0] == '#' ) { in.getline(tmpLine,1000); }
	stringstream s;
	s << tmpLine;
	s >> width >> height;
	size = width * height;
	tmpLine[0] = '#';
	while( tmpLine[0] == '#' ) { in.getline(tmpLine,1000); }
	if (strcmp(tmpLine, "255") != 0)
	  msgs::error("Can't understand the header");
	image_string = new unsigned char[size];
	in.read((char*)image_string,size);
}

void GRImage::saveTo(const string& file_name)
{
  ofstream out(file_name.c_str(), ios::out | ios::binary);

  int pl = file_name.length() - 1;
	char ext[4];
	int count = pl;
	while (file_name[count] != '.'){
		ext[pl - count] = file_name[count];
		count--;
	}
	ext[pl - count] = '\0';

	if (strcmp(ext,"wb") == 0 )
		printOnBW(out);
	else if(strcmp(ext,"mgp") == 0 )
		printOnPGM(out);
	out.close();	
}


// -------------------- LUTImage -----------------------//

int GRLUTImage::getPixel(int i, int j) const
{
	return   theLookUpTable.get(GRImage::getPixel(i,j));
}

int GRLUTImage::getPixel(int n) const
{
  	return   theLookUpTable.get(image_string[n]);
}


void GRLUTImage::printOn( ostream& out)const
  // Output the image into a stream
{
  
  int w = width;
  int h = height;
  unsigned char p[8];
   
  p[0] = ((w>>24) & 0xff);	/* save 4-byte width  */
  p[1] = ((w>>16) & 0xff);
  p[2] = ((w>> 8) & 0xff);
  p[3] = ( w      & 0xff);
  
  p[4] = ((h>>24) & 0xff);	/* save 4-byte height */
  p[5] = ((h>>16) & 0xff);
  p[6] = ((h>> 8) & 0xff);
  p[7] = ( h      & 0xff);

  out << p;

  for (int i=0;i<size;i++)
    out << theLookUpTable.get(image_string[i]);
}



////////////////////////////////////////////////
//
//  CIMAGE
//
///////////////////////////////////////////////


void convert(const CImage* ci, GRImage* gi)
{
  for (int i=0;i<ci->getSize();i++){
    int Y = int(0.299*ci->getRedPixel(i) + 0.587*ci->getGreenPixel(i)+0.114*ci->getBluePixel(i));
    gi->setPixelCliped(i, Y );
  }
}

CImage::CImage( const string& in_file_name, FILE_TYPE ft):AImage( in_file_name )
  // Create image from the input file
{
  ifstream in(in_file_name.c_str(),ios::in |ios::binary);
  if (ft == PPM )
    readFromPPM(in);
  in.close();
}

void CImage::printOnPPM( ostream& out)const
// Output the image into a stream
{
  out << "P6" << endl << "#MAD Alex" << endl<< width <<" "<<height<<endl <<255 << endl;
  
  out.flush();
  for (int i=0;i<size;i++){
    //	out << getRedPixel(i) << " " << getGreenPixel(i) << " " << getBluePixel(i) << endl;
    unsigned char p[3];
    p[0] = ( unsigned char)getRedPixel(i);
    p[1] = ( unsigned char)getGreenPixel(i);
    p[2] = ( unsigned char)getBluePixel(i);
    out.write((char*)p,3);
  }
  out.flush(); 
}

void CImage::readFromPPM(istream &in)
{
  char tmpLine[1000];
  in.getline(tmpLine,1000);
  if ((strcmp(tmpLine, "P3") == 0)){
    tmpLine[0] = '#';
    while( tmpLine[0] == '#' ) { in.getline(tmpLine,1000); }
    stringstream s;
    s << tmpLine;
    s >> width >> height;
    size = width * height;
    greenImage = blueImage = redImage = GRImage( width, height );
    
    tmpLine[0] = '#';
    while( tmpLine[0] == '#' ) { in.getline(tmpLine,1000); }
    //error_check(strcmp(tmpLine, "255") != 0, "Can't understand the header");
    
    for (int i=0;i<size;i++){
      int r,b,g;
      in >> r >> g >> b;
	    setRedPixel(i,r); setGreenPixel(i,g); setBluePixel(i,b);
    }		
  }else if ((strcmp(tmpLine, "P6") == 0)){
    tmpLine[0] = '#';
    while( tmpLine[0] == '#' ) { in.getline(tmpLine,1000); }
    stringstream s;
    s << tmpLine;
    s >> width >> height;
    size = width * height;
    greenImage = blueImage = redImage = GRImage( width, height );
    
    tmpLine[0] = '#';
    while( tmpLine[0] == '#' ) { in.getline(tmpLine,1000); }
    //error_check(strcmp(tmpLine, "255") != 0, "Can't understand the header");
    
	  for (int i=0;i<size;i++){
	    unsigned char p[3];
	    in.read((char*)p,3);
	    setRedPixel(i,p[0]); setGreenPixel(i,p[1]); setBluePixel(i,p[2]);
	  }	
  } else if (strcmp(tmpLine, "P5") == 0){
    tmpLine[0] = '#';
    while( tmpLine[0] == '#' ) { in.getline(tmpLine,1000); }
    stringstream s;
    s << tmpLine;
    s >> width >> height;
    size = width * height;
    greenImage = blueImage = redImage = GRImage( width, height );
    
    tmpLine[0] = '#';
    while( tmpLine[0] == '#' ) { in.getline(tmpLine,1000); }
    //error_check(strcmp(tmpLine, "255") != 0, "Can't understand the header");
    
    for (int i=0;i<size;i++){
      unsigned char p[1];
      in.read((char*)p,1);
      setRedPixel(i,p[0]); setGreenPixel(i,p[0]); setBluePixel(i,p[0]);
    }
  }else
    cout << "Couldn't read an input image";
}

void CImage::saveTo(const string& file_name)
{
  //	ofstream out(file_name, ios::out | ios::binary, filebuf::sh_write);
  ofstream out(file_name.c_str(), ios::out | ios::binary);
  
  if( getFileType(file_name) == PPM )
    printOnPPM(out);
  else
    msgs::error("Cannot recognise image file format. Please add an extension to the file name. Supported extensions: .ppm ");
  out.close();	
}

