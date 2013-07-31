#ifndef PROPOSEDTRANSCRIPT_H
#define PROPOSEDTRANSCRIPT_H

#include "Centers.h"

class ProposedTranscript: public Transcript{
 public:
  void to_z(vector<Center> &v, unsigned int &size, Arguments &A);
  ProposedTranscript(unsigned int &my_t, Centers* centers, Arguments &A);
};

#endif
