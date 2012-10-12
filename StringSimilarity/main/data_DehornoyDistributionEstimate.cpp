
// Contents: Generates pairs of words and computes similarity measures of a given type.
//           Estimates the probasbility distribution of the values for these measures
//           and records it in the file as an array of values (no bin counts for now)
//           
// Principal Author: Alexei Miasnikov (2004)
//
// Status:
//


#include <values.h>
#include "global.h"
#include "FreeGroup.h"
#include <sstream>
#include "RanlibCPP.h"
#include "RandomWord.h"
#include <strstream>
#include "SimilarityMeasures.h"
#include <vector>
#include <algorithm>
#include "Map.h"
#include "FRDFIT.h"
#include "ConfigFile.h"

int main(int argc, char** argv)
{
  Map m;
  
  cout << "## Estimate a distribution of a measure between a random word and its Dehornoy form. Alexei Miasnikov (2004)" << endl << endl; 
  if (argc != 2){
    cout << "Usage: executable <configuration>" << endl;
    exit(0);
  }
  
  //  strstream s;
  //  s << argv[1] << endl;
  //  FreeGroup F;
  //  s >> F;  
  

  ConfigFile cfg(argv[1]);
  cout << "Configuration : " << endl
       << cfg << endl << endl;
  
  
  int nGens =  cfg.getValue("NGENS");
  int number = cfg.getValue("SAMPLE_SIZE");
  int wordLen = cfg.getValue("WORDS_LENGTH");
  string measure_type = string(cfg.getValue("MEASURE"));
  string scrambler = string(cfg.getValue("SCRAMBLER"));


  FreeGroup F(nGens);


  UniformRandomWord urw( nGens,wordLen  -wordLen / 2,wordLen + wordLen/2 );

  // Initialize the similarity measure
  StringSimilarityMeasure* sm;
  if  (measure_type == string("HAMMING")) sm = new HammingDistance();
  else if (measure_type == string("CYCLIC_HAMMING")) sm = new SubwordHammingDistanceCyclic();
  else if (measure_type == string("WHITEHEAD_GRAPH")) sm = new WhiteheadGraphSimilarity(F); 
  else if (measure_type == string("EDITING")) sm = new EditingDistance();
  else if (measure_type == string("CYCLIC_EDITING")) sm = new SubwordEditingDistanceCyclic();
  else if (measure_type == string("MOTIVE")) {
    if (nGens > 26){
      cout << "TMP hack. Do not allow more then 26 gens" << endl;
      exit(0);
    }
    
    VectorOf<Chars> gNames(nGens);
    for (char g=0;g<nGens;g++){
      gNames[g] = Chars(char('A'+g) );
    }
    FreeGroup nF( gNames );
    
    cout << nF << endl;
    sm = new MotiveSimilarity(nF);
    
  }
  else {
    msgs::error("Cannot recognize measure type.");
  }
  
  
  // Vector of measure values
  vector<double> measureDistr(number);
  
  // generate "number" of pairs
  for ( int lCount=0;lCount<number;lCount++ ){
    
    // Randomly generate a pair of a wortd and a Dehornoy form of another word
    Word w1 = urw.word( );
    Word w2 = urw.word( );

    /*
      FRDFIT frdfit(F);
      int nh = 0;
      Word dw = frdfit.compute(w2,nh);
    */
    
    DehornoyForm df( nGens , w2 );
    Word dw = df.getDehornoyForm( );
    
    measureDistr[lCount] = sm->measure( w1,dw );
  }
  
  sort(measureDistr.begin(),measureDistr.end());
  for ( int lCount=0;lCount<number;lCount++ ){
    cout << measureDistr[lCount] << endl;
  }
  
  int n_repeat = 100;
  
  // TMP Do  evaluation on elements with their Dehornoy forms
  cout << "Evaluate on TRUE DEHORNOY FORMS --------------------------- " << endl;
  double sum_pvalue = 0.0;
  double sum_pvalue_squared = 0.0;

  for (int i=0;i<n_repeat;i++){
    Word w1 = urw.word( wordLen );
    FRDFIT frdfit(F);
    int nh = 0;
    Word dw1 = frdfit.compute(w1,nh);
    
    double conf_value =  getLeftTaleConfidenceValue(measureDistr,sm->measure( w1,dw1 ));
    cout << conf_value << endl ;
    sum_pvalue += conf_value;
    sum_pvalue_squared += conf_value*conf_value;
    
  }
  
  cout << endl << "Accuracy TRUE DEHORNOY FORMS : " << sum_pvalue  / double(n_repeat) 
       << " +/- " << sqrt((sum_pvalue_squared - sum_pvalue*sum_pvalue / double(n_repeat)) / double(n_repeat - 1))<< endl;


  // TMP Do  evaluation on elements with independent Dehornoy forms
  cout << "Evaluate on INDEPENDENT DEHORNOY FORMS --------------------------- " << endl;
  sum_pvalue = 0.0;
  sum_pvalue_squared = 0.0;

  
  for (int i=0;i<n_repeat;i++){
    Word w1 = urw.word( wordLen );
    Word w2 = urw.word( wordLen );
    FRDFIT frdfit(F);
    int nh = 0;
    Word dw2 = frdfit.compute(w2,nh);
   
    double conf_value = getLeftTaleConfidenceValue(measureDistr,sm->measure( w1,dw2 ));
    cout << conf_value << endl ;
    sum_pvalue += conf_value;
    sum_pvalue_squared += conf_value*conf_value;

  }
  
  cout << endl << "Accuracy INDEPENDENT DEHORNOY FORMS : " << sum_pvalue / double(n_repeat) 
       <<  " +/- " << sqrt((sum_pvalue_squared - sum_pvalue*sum_pvalue / double(n_repeat)) / double(n_repeat - 1)) << endl;
  return 0;
}
