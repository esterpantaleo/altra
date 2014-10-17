#ifndef SUMMARY_H
#define SUMMARY_H

#include "Arguments.h"

//if transcripts are +,100110,110100,101010 -,000001,010000
//then the transcript model z is 100110110100101010000001010000
//i.e., the concatenated bitsets
//this maps a transcript model z to aa loglikelihood vector
typedef std::map<boost::dynamic_bitset<>,vector<double> > SummaryAsMap; 

//this is a pair with the average loglikelihood and a pair (transcript model z and number of states over which the average loglikelihood was computed
typedef std::pair<double,std::pair<boost::dynamic_bitset<>, unsigned int> > SummaryAsPair;

class Summary{
public:
  //the total number of steps of the MCMC (after burnin)
  unsigned int number_steps;
  vector<SummaryAsPair> myPairs;

  Summary(SummaryAsMap * myMap);
  //this function is used by print
  void print_util(std::vector<SummaryAsPair>::iterator it, ofstream & ssfile, Arguments &A);
  //this function prints the Summary (a list of transcript models
  //ordered by the loglikelihood value) and the GenePred file
  void print(Arguments &A);
};

#endif 
