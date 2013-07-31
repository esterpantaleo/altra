#include "ProposedTranscript.h"

void ProposedTranscript::to_z(vector<Center> &v, unsigned int &size, Arguments &A){
  //z must be all zeros                                                
  for (unsigned int c = 0; c < size; c ++){
    if (v[c].is_used_proposed[t])//==1               
      for (unsigned int cc = (v[c].v3_proposed[t] -> c); cc < (v[c].v5_proposed[t] -> c); cc ++)
	z[A.lC - cc - 2] = 1;
  }
}

ProposedTranscript::ProposedTranscript(unsigned int &my_t, Centers* centers, Arguments &A):
  Transcript(my_t, A){
  to_z(centers -> v, centers -> size, A);
}
