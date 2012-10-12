/*
 *   $Id: MotivePatternWrapper.cpp,v 1.1 2005/11/28 23:25:56 amiasnik Exp $
 */
 
// Contents: 
//
// Principal Author: Alexei Miasnikov (2004)
//
// Status:
//
// Revision History:
//


const char upperCaseShift  = 'a' - 'A';

#include "global.h"
#include "Word.h"
#include "MotivePatternWrapper.h"
#include <values.h>
#include "Chars.h"
#include <sstream>

const char* mpe_executable = "./mpe_exec";

void MotivePatternWrapper::initializeSequence( const vector<Word>&  input)
{
  input_file = string(tmpnam(NULL));
  //cout << "Input file name : " << input_file << endl;
  ofstream input_stream( input_file.c_str() );
  for (int i=0;i<input.size();i++){
    input_stream << ">W" << i+1 << "  " << i+1 << endl;
    for (int j=0;j<input[i].length();j++){
      Chars gen = theGroup.nameOfGenerator( abs(input[i][j])-1 );
      if ( ord(input[i][j]) > 0 )
	input_stream  <<  gen;
      else
	for (int n=0;n<gen.length();n++){
	  if (gen[n] >='a' && gen[n] <= 'z')
	    input_stream  << char(gen[n] - upperCaseShift);
	  else
	    input_stream  << char(gen[n] + upperCaseShift);
	}
	
    }
  
    input_stream << endl;
  }
  input_stream.close();
}

bool MotivePatternWrapper::compute()
{
  output_file = string(tmpnam(NULL));

  stringstream s;
  s << mpe_executable << " -l2 -w2 -c1 -k2 -v -i" << input_file << " -o" << output_file << " > NULL " <<  '\0' << flush;

  //cout << s.str().c_str() << endl;
  
  int ret = system(s.str().c_str());
  //  s.free();
  if ( ret < 0 ){
    cout << "Cannot execute motive pattern extractor..." << endl;
    exit(0);
  }
  
}

vector<MotivePattern> MotivePatternWrapper::getPatterns()
{
  vector<MotivePattern> patterns;

  MotivePattern mp;
  ifstream pr_stream(output_file.c_str());

  while (!pr_stream.eof()){
    
    mp.occurrences = 0;
    mp.sequences  = 0;
    string pattern;
    
    string tmp_s;
    pr_stream >> tmp_s;
    if ( (mp.occurrences = atoi(tmp_s.c_str())) ){
      pr_stream >> tmp_s;
      mp.sequences = atoi(tmp_s.c_str());
      
      pr_stream >> pattern;

      stringstream s;
      for (int i=0;i<pattern.length();i++)
	s << pattern[i] << " ";
      
      Chars errMsg;
      Word w = theGroup.readWord(s,errMsg);
      if (errMsg.length() > 0){
	cout << "Couldn't read a pattern from file" << endl;
	exit(0);
      }
      
      
      mp.thePattern = w;
      
      //    cout << occurances << " " << sequences << " ";
      //      theGroup.printWord(cout,w);
      //      cout << endl;
      patterns.push_back( mp );
    }
    
  }
  pr_stream.close();
  //cout << "Patterns red successfully" << endl;
  
  return patterns;
}



