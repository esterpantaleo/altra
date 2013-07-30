#ifndef TRANSCRIPT_H
#define TRANSCRIPT_H

#include "Center.h"

class Transcript{
public:
  unsigned int t;
  boost::dynamic_bitset<> z;
  boost::dynamic_bitset<> modifiedExon;

  bool operator< (const Transcript &rhs) const;
  bool operator== (const Transcript &rhs) const;
  virtual void to_z(vector<Center> &v, unsigned int &size, Arguments &A);

  Transcript(){}  
  Transcript(unsigned int &my_t, Arguments &A);
  Transcript(Arguments &A);
  Transcript(unsigned int &my_t, vector<Center> &v, unsigned int &size, Arguments &A);
  int len(Arguments &A);
};

#endif
