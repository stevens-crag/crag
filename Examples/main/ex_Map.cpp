// Copyright (C) 2006 Aleksey Myasnikov
// Contents: Example for class Map
//
// Principal Authors:  Aleksey Myasnikov
//
// Revision History:
//

#include "Map.h"
#include "Alphabet.h"

int main( )
{

  // Read domain's alphabet from the inout stream
  FiniteAlphabet domain;
  cout << "Input alphabet of the domain ( {a,b,c, ...} ): " << flush;
  cin >> domain;

  // Input range
  FiniteAlphabet range;
  cout << "Input alphabet of the range ( {x,y,z, ...} ): " << flush;
  cin >> range;  
  
  // Create a set of random images 
  vector<Word> random_images(domain.size());
  for (int i=0;i<random_images.size();i++)
    random_images[i] = Word::randomWord(range.size(),2,3);


  //&Maps; How do I create a map 
  
  // create a map from domain to range, defined by random_images
  Map m(domain,range, random_images);
  
  
  
  //& Maps; How do I print a map
  
  // Print a map in the standard output
  cout << "Map from " << domain << " to " << range << ":" << endl << m << endl; 
  


  //& Maps; How do I input a map

  // Read a map from the standard input stream
  
  cout << "Enter a map from " << domain << " into " << range << endl
       << "For example if domain={a,b,c}  and range={x,y,z} : { a -> x y , b -> y z , c -> x z }" << endl;
  cin >> m;
  cout << "You entered : " << m << endl;
  

  //&Maps; How do I compute a word's image

  cout << "Enter a word in " << domain << endl;
  Word w = domain.readWord( cin );  

  // compute the image
  Word w_image = m.imageOf( w );
  
  // Print results
  cout << "Map results : ";
  domain.printWord(cout,w);  cout << " -> ";  
  range.printWord(cout,w_image); cout << endl;
  

  //&Maps; How do I compute a composition of two maps
  Map m1( domain,range );
  Map m2( range,domain );
  
  cout << "Enter a map from " << domain << " into " << range << ":" << flush;
  cin >> m1;

  cout << "Enter a map from " << range << " into " << domain << ":" << flush;
  cin >> m2;

  cout << "The composition " << composition( m1,m2 ) << endl;

  

  //&Maps; How do I generate a random automorphism
  
  // Generate a random automorphism by composing 10 Whitehead automorphisms
  // in a free group of rank 5
  cout << "Random auto in F_5 of length 3 " << RMap::getRandomAuto( 5,3 ) << endl;

  //&Maps; How do I generate a random Whitehead's automorphism

  cout << "Random Whitehead auto in F_5 : " << RMap::getRandomWhiteheadAuto( 5 ) << endl;
  
  return 0;
}
