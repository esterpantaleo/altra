#ifndef CENTERS_H
#define CENTERS_H

#include "Transcript.h"

class Centers{
public:
  unsigned int size, K;
  string strand;
  vector<SS> three_primes, five_primes;
  vector<Center> v;// contains the set of centers, of size "size"          
  vector<Transcript> transcripts;
  double ar_t, ar_den_t;
  vector<double> ar, ar_den; 

  void do_nothing(){}
  void clear();
  void print_ar(string strand, Arguments &A);
  void clear_ar(Arguments &A);
  void init_proposal();
  void init(string my_strand, vector<unsigned int> &pos3, vector<unsigned int> &pos5, Arguments &A);
  bool get_c3_c5(unsigned int &genePred5SS, unsigned int &genePred3SS, unsigned int &c3, unsigned int &c5, Arguments &A);
  void init2(vector<unsigned int> &genePred3SS, vector<unsigned int> &genePred5SS, Arguments &A);
  bool check_init(string &genePredFile, Arguments &A);
  void init(string &genePredFile, Arguments &A);
};

#endif
