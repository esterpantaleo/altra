#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include "utils.h"

class Arguments{
 public:
  //lC=number of elements of listCoordinate                       
  unsigned int RL, K, powK, pK, nK, N, MC_STEPS, MC_BURNIN, MC_EQ, MC_THIN, lC, genePrediction, MIN_EX_LEN, MAX_EX_LEN_DIV_2, MIN_IN_LEN, MAX_IN_LEN;
  unsigned long int rSeed;
  //LC=length of the region, Csum=sum of the C, parameters of the priors on log_lambda and epsilon   
  double LC, sumC, alphaSS, Csum, mu_e, var_e, sd_e, sd_q, var_l, ar_e, ar_l, scale_junction_count_by, NORMALIZATION, a, b, c;
  string genePredFile, chr, genePredOut, utils, output2summary, run_output2summary, outfilellk, outfile1, argv3, argv14;
  vector<string> my_types;
  ostringstream buffer;
  ofstream outdatallk, outdata1;
  vector<unsigned int> L, pos3, pos5, neg3, neg5;
  vector<double> C, log_C, listCoordinate;
  gsl_matrix* Z;//a matrix with rows the binary form of z (ulong to dynamic_bitset from z=0 to z=A.powK ) 
  gsl_matrix* VAR_q;
  gsl_matrix* VAR_l;
  //random seed                                 
  gsl_rng* rnd;
  const gsl_rng_type* Tor;

  Arguments(int argc,const char* argv[]);
  ~Arguments();

  void set0(unsigned int &my_pK, unsigned int &my_nK);
  void set(unsigned int &my_pK, unsigned int &my_nK);
  void set2(unsigned int &my_pK, unsigned int &my_nK);
  
  void clear_ar();
  void print_ar(const unsigned int &num_steps);
  void print_ar_my_update(const unsigned int &num_steps);
};
    
#endif //ARGUMENTS_H 
