
#ifndef _MAP_H_
#define _MAP_H_

#include <iostream>

#include <vector>
#include "Word.h"
#include "errormsgs.h"

////////////////////////////////////////////////////////////////////
//                                                                //
// Map                                                            //
//                                                                //
////////////////////////////////////////////////////////////////////


//! Defines a representation of a map between two sets of words.
/*!
  Both the domain and the range of a map are defined by finite alphabets. 
 */
class Map
{
  
public:
  
  //! Default constructor.
//   Map( ) : 
//     theDomain(0), 
//     theRange(0),
//     theDomainAlphabet(&(InfiniteAlphabet::defaultAlphabet)),
//     theRangeAlphabet(&(InfiniteAlphabet::defaultAlphabet)) {} 

  // copy constructor supplied by compiler
  

  
  //! Constructor. Creates a trivial map.
  /*! If no  images of letters of the domain are
    specified the map, by  default, becomes the trivial map, i.e.
    the map mapping everything to the identity.
    \param domain_gens - the size of the domain's alphabet
    \param range_gens - the size of the range's alphabet
  */
  Map( int domain_gens, int range_gens ) 
    : theDomain(domain_gens),
    theRange(range_gens),
    theDomainAlphabet(domain_gens),
    theRangeAlphabet(range_gens),
    theGeneratingImages(domain_gens)  { }

  
  //! Constructor. Creates a map defined by the  images of an alphabet.
  /*! Constructs a map, given a domain and range group, and images
    for the generators of the domain
    \param domain_gens - the size of the domain's alphabet
    \param range_gens - the size of the range's alphabet
    \param generatingImages - vector of images of letters of the domain.
  */
  Map(  int domain_gens, int range_gens, const vector<Word>& generatingImages ) 
    :theDomain(domain_gens),
    theRange(range_gens),
    theDomainAlphabet(domain_gens),
    theRangeAlphabet(range_gens),
    theGeneratingImages(generatingImages)  
  { 
    
    if (generatingImages.size() != domain_gens)
      msgs::error("wrong number of generating images in  Map::Map(domain,range,generatingImages)");
    
  }

  //! Constructor. Creates a trivial map.
  /*! If no  images of letters of the domain are
    specified the map, by  default, becomes the trivial map, i.e.
    the map mapping everything to the identity.
    \param domain_alph - alphabet of the domain
    \param range_alph - alphabet of  the range
  */
  Map( const FiniteAlphabet&  domain_alph, const FiniteAlphabet& range_alph ) 
    : theDomain(domain_alph.size()),
      theRange(range_alph.size()),
      theDomainAlphabet(domain_alph),
      theRangeAlphabet(range_alph),
      theGeneratingImages(domain_alph.size())  { }

  
  //! Constructor. Creates a map defined by the  images of an alphabet.
  /*! Constructs a map, given a domain and range group, and images
    for the generators of the domain
    \param domain_alph - alphabet of the domain
    \param range_alph - alphabet of  the range
    \param generatingImages - vector of images of letters of the domain.
  */
  Map(  const FiniteAlphabet&  domain_alph, const FiniteAlphabet& range_alph, 
	const vector<Word>& generatingImages ) 
    :theDomain(domain_alph.size()),
     theRange(range_alph.size()),
     theDomainAlphabet(domain_alph),
     theRangeAlphabet(range_alph),
     theGeneratingImages(generatingImages)  
  { 
    
    if (generatingImages.size() != theDomain)
      msgs::error("wrong number of generating images in  Map::Map(domain,range,generatingImages)");
    
  }
  // destructor supplied by compiler

  ////////////////////////////////////////////////////////////////////
  //                                                                //
  // Standard Operators                                             //
  //                                                                //
  ////////////////////////////////////////////////////////////////////

  //! Comparison operator.
  /*!
    Two maps are equal if the sizes of domains and ranges coincide respectively and
    images are defined by the same words.
   */
  bool operator == ( const Map& m ) const 
  {
    return ( theDomain == m.theDomain && theRange == m.theRange && theGeneratingImages == m.theGeneratingImages );
  }

  ////////////////////////////////////////////////////////////////////
  //                                                                //
  // Accessors / Modifiers                                          //
  //                                                                //
  ////////////////////////////////////////////////////////////////////

  //! Returns the size of the domain's alphabet
  int domainSize( ) const { return theDomain; }

  //! Returns  domain's alphabet
  const FiniteAlphabet&  domainAlphabet( ) const { return theDomainAlphabet; }

  //! Returns the size of the range's alphabet
  int rangeSize( ) const { return theRange; }

  //! Returns  range's alphabet
  const FiniteAlphabet&  rangeAlphabet( ) const { return theRangeAlphabet; }

 
  //! Returns the vector of images 
  const vector<Word>& generatingImages( ) const { return theGeneratingImages; }

  //! Return the  image of the ith letter
  Word generatingImages( int i ) const { 
    if (i < 0 || i >= theDomain)
      msgs::error("generating image index out of range in  Map::generatingImages(int)");
    return theGeneratingImages[i]; 
  } 

 
  //! Sets or modifies generating images
 void setGeneratingImages( const vector<Word> gi ) 
  {
    if (gi.size() != theDomain)
      msgs::error("wrong number of generating images in  Map::setGeneratingImages( const VectorOf<Word> )");
    theGeneratingImages = gi;
  }

  
 //!  Assigns to the i-th (0-based) generating image 
 /*!
   \param i - the  index of a letter
   \param e - the new image
 */
 void setGeneratingImages( int i, const Word& e  ) 
  {

    if (i < 0 || i >= theDomain)
      msgs::error("generating image index out of range in  Map::setGeneratingImages(int, const Word&)");
    theGeneratingImages[i] = e;
  }


	
  ////////////////////////////////////////////////////////////////////
  //                                                                //
  // Mapping methods                                                //
  //                                                                //
  ////////////////////////////////////////////////////////////////////

  // computing images:
  // formally:

 //! Compute the image of a word
 /*!
   Takes a formal word in the generators of domain and evaluates its
   `image' in range. The image is computed as follows. Let input word
   \f$ w = x_1 x_2 \ldots x_n,\f$ where \f$x_i \in \f$ Domain and 
   \f$image(x) \in \f$ Range be the image of the letter \c x, then
   image of \c w - \f$image(w) = image(x_1) image(x_2) \ldots image(x_n)\f$
   \param w - a word in the domain
   \return an image of \c w.
   
 */
  Word imageOf( const Word& w ) const;


  //! Map-theoretic composition
  /*!
    Returns the composition of two maps
    \param firstMap - the first map
    \param secondMap - the second map
    \return return the map constructed as a composition \c secondMap(\c firstMap)
   */
  friend Map composition( const Map& firstMap,
			  const Map& secondMap );
   

  //! Map composition operator
  /*!
    See \link composition() composition( ... ) \endlink for details.
   */
  Map operator | ( const Map& secondMap )
  { return composition(*this, secondMap); }
  


  ////////////////////////////////////////////////////////////////////
  //                                                                //
  // I/O                                                            //
  //                                                                //
  ////////////////////////////////////////////////////////////////////


  //------------------ Related global functions --------------------//
  
  //! Stream output operator
  friend ostream& operator << ( ostream& out, const Map& m ){
    m.printOn(out);
    return out;
  }

  //! Stream input operator
  friend istream& operator >> ( istream& in, Map& m ){
    m.readFrom(in);
    return in;
  }

  //! Stream output
  void printOn(ostream& ostr) const;
 
  //! Stream input
  void readFrom(istream& in);
private:
  char readChar(istream& in)const;


  // real data:
  int theDomain;
  int theRange;
  vector<Word> theGeneratingImages;

  FiniteAlphabet theDomainAlphabet;
  FiniteAlphabet theRangeAlphabet;
  
};


////////////////////////////////////////////////////////////////////
//                                                                //
// Random Maps                                                    //
//                                                                //
////////////////////////////////////////////////////////////////////
  
//! Namespace for random map generating methods
namespace RMap 
{
  //! Returns a randomly generated Whitehead automorphism
  /*!
    \param n - the rank of a free group
    \return returns a randomly generated Whitehead automorphism
  */
  Map getRandomWhiteheadAuto(  int n );

  //! Returns a randomly generated automorhism of a given length
  /*!
    \param n - the rank of a free group
    \param l - the number of Whitehead automorphisms 
   */
  Map getRandomAuto( int n, int l);
};

#endif
