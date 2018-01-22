
#include "Map.h"
#include "RanlibCPP.h"

//------------------------------------------------------------------//
//----------------------------- Map --------------------------------//
//------------------------------------------------------------------//


Word Map::imageOf( const Word& w ) const {
  Word image;
  
  for ( auto I = w.begin(); I!=w.end(); I++){
    auto g = *I;
    image *= ( g > 0 ? theGeneratingImages[abs(g)-1] :  theGeneratingImages[abs(g)-1].inverse() );
  }

  return image; 
  
}

void Map::printOn(ostream& ostr) const
{
  ostr << "{ ";
  int len = generatingImages().size();
  for( int i = 0; i < len; ++i ) {	 
    if ( i > 0 ) ostr << ", ";
    theDomainAlphabet.printWord(ostr,Word( i+1));
    ostr << " -> ";
    theRangeAlphabet.printWord(ostr, generatingImages(i));
  }
  ostr << " }";
}

void Map::readFrom(istream& in) 
{

  // clear the images
  for (int i=0;i<theDomain;i++)
    theGeneratingImages[i] = Word();

  // read the first symbol which supposed to be '}'
  char ch = readChar( in );
  
  if (ch!='{') 
    msgs::error("Map::readFrom(): syntax error.");
  
  int rangeCount = 0;
  while (1){
    
    if ( ++rangeCount > theRange )
      msgs::error("Map::readFrom(): Range dimension exceeded. Syntax error.");

    //cout << rangeCount << ": "  <<endl;
    Parser pg(in,&theDomainAlphabet);
    pg.parse();
    //cout << "Generator parsed" << endl;
    Word g(pg.getWord());
    //cout << "Gemerator parsed" << endl;

    //cout << "Read generator"  <<endl;

    if ( g.length() != 1 )
      msgs::error("Map::readFrom(): Not a generator. Syntax error.");
    
    if (pg.getWordTerminalSymbol() != '-' )
      msgs::error("Map::readFrom(): syntax error.");

    //    delete pg;

    ch = readChar( in );
    if (ch!='>') 
      msgs::error("Map::readFrom(): syntax error.");

    Parser pi(in,&theRangeAlphabet);
    pi.parse();
    Word image_g(pi.getWord());

    //cout << "Read image"  <<endl;
  
    auto gI = g.begin();

    int g_i  = *(gI)-1;

    //cout << "Iterator " << *gI << " - " << g_i << endl;

    
    theGeneratingImages[g_i] = image_g;
    
    if (pi.getWordTerminalSymbol() == '}'  ) break;
    if (pi.getWordTerminalSymbol() != ',' && pi.getWordTerminalSymbol() != ';'){
      msgs::error("Alphabet::readVector(): syntax error.");
    }
    //    delete pi;
  }
  
}	

char Map::readChar( istream& in)const{
  char ch = ' ';
  while ( ch == ' ' || ch=='\t' || ch == '\n' || ch == '\r' ){
    in >> ch;
  } 
  return ch;
}

Map RMap::getRandomWhiteheadAuto(  int nOfGens )
{
  vector<Word> image(nOfGens);
  int ai = RandLib::ur.irand(0,nOfGens-1);
  auto a = ai+1;
  if (RandLib::ur.rand() < 0.5)
    a = -a;
  image[ai] = Word(a);
  for (int i=0;i<nOfGens;i++){
    if (i != ai){
      int choice = RandLib::ur.irand(0,4);
      if (choice == 0)
	image[i] = Word(i+1);
      else if (choice == 1)
	image[i] = Word(-(i+1));
      else if (choice == 2)
	image[i] = Word(i+1)*Word(a);
      else if (choice == 3)
	image[i] = Word(-a)*Word(i+1);
      else
	image[i] = Word(-a)*Word(i+1)*Word(a);
    }
  }

  Map m(nOfGens,nOfGens,image);
  return m;
  
}

Map RMap::getRandomAuto( int nOfGens, int len )
{
  Map m = getRandomWhiteheadAuto(nOfGens);
  for (int i=1;i<len;i++)
    m = composition( m, getRandomWhiteheadAuto( nOfGens ));

  return m;
}


//---------------------------- global -----------------------------//


Map composition( const Map& firstMap, const Map& secondMap )
{
  
  int i;
  
  if (secondMap.domainSize() != firstMap.rangeSize())
    msgs::error("tried to compose differing domain and range in"
	  " ::composition(const Map&, const Map&)");
  
  vector<Word> images(firstMap.generatingImages().size());
  for (i=0; i < images.size(); i++) {
    Word w = firstMap.generatingImages(i);
    images[i] = secondMap.imageOf( w );
  }
  
  return Map(firstMap.domainAlphabet(),secondMap.rangeAlphabet(),images);
}


