#ifndef CENTER_H
#define CENTER_H

#include "Arguments.h"

class SS{
public:
  unsigned int c, n, n_proposed;//c = the element in the list of coordinates (pos3 or pos5); n = number of transcripts that use that splice site 
 
  SS(unsigned int &cc);
};

bool SS_less(const SS* lhs, const SS* rhs);
bool SS_equal(const SS* lhs, const SS* rhs);

class SSpointers{
public:
  vector<SS*> v;
  unsigned int size;
};

class Center{
public:
  SSpointers ss3Pointers, ss5Pointers;
  vector<SS*> v5, v3, v5_proposed, v3_proposed;//a vector of K-dim: my_3 is the list of 3' ss used by each of the K transcripts at center c             
  unsigned int n, n_proposed;// number of transcripts that use center, total number of transcripts  
  vector<bool> is_used, is_used_proposed; // vector of size K is_used[k]==1(0) means center is used by transcript

  void clear();
  Center(unsigned int &i3, unsigned int &j5, vector<SS> &ss3, vector<SS> &ss5, Arguments &A);
  bool get(unsigned int &c3, unsigned int &c5, SS*& s3, SS*& s5);
  SS* draw_SS_from_posterior(SSpointers &ssPointers, Arguments &A);
  void init_not_used(Arguments &A);
};

class ModifiedCenter{
public:
  unsigned int t, c;

  ModifiedCenter(){}
  ModifiedCenter(unsigned int &my_c, unsigned int &my_t);
};

bool is_less_than(ModifiedCenter x, ModifiedCenter y);
bool is_equal_to(ModifiedCenter x, ModifiedCenter y);

#endif
