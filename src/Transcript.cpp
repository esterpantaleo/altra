#include "Transcript.h"

bool Transcript::operator< (const Transcript &rhs) const{
  return (z < rhs.z);
}

bool Transcript::operator== (const Transcript &rhs) const{
  return (z == rhs.z);
}
  
void Transcript::to_z(vector<Center> &v, unsigned int &size, Arguments &A){
  for (unsigned int c = 0; c < size; c ++){
    if (v[c].is_used[t]){//==1
      for (unsigned int cc = (v[c].v3[t] -> c); cc < (v[c].v5[t] -> c); cc ++)
	z[A.lC - cc - 2] = 1;
    }
  }
}

Transcript::Transcript(unsigned int &my_t, Arguments &A){
  t = my_t;
  for (unsigned int i = 0; i < A.lC - 1; i++)
    z.push_back(0); 
}

Transcript::Transcript(Arguments &A){
  for (unsigned int i = 0; i < A.lC - 1; i++)
    z.push_back(0);                                      
}

Transcript::Transcript(unsigned int &my_t, vector<Center> &v, unsigned int &size, Arguments &A){
  t = my_t;
  for (unsigned int i = 0; i < A.lC - 1; i++)
    z.push_back(0);                       
  to_z(v, size, A);
}

int Transcript::len(Arguments &A){
  int Len = 0;
  
  for (unsigned int c = 0; c < A.lC - 1; c ++){
    if (z[A.lC - c - 2])//==1                                  
      Len += A.L[c];
  }
  if (Len == 0) 
    return 0;
  else 
    Len -= A.RL - 1;
  //Len = max(0, Len);
  
  return Len;
}    

