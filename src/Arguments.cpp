/** \file Arguments.cpp
 *

 * Copyright (C) 2012-2013 Ester Pantaleo


 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by


 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of


 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.


 */

#include "Arguments.h"

Arguments::~Arguments(){
  if (VERBOSE==1)
    outdatallk.close();
  gsl_matrix_free(Z);
  gsl_matrix_free(VAR_q);
  gsl_matrix_free(VAR_l);
  gsl_rng_free(rnd);
}

void Arguments::clear_ar(){
  ar_l = 0.;
  ar_e = 0.;
}

void Arguments::print_ar(const unsigned int &num_steps){
  cout<<"ACCEPTANCE RATIO FOR LAMBDA=                    " << ar_l << "/" << N * num_steps  << "=" << ar_l / (double) (N * num_steps) << "\n";
  cout <<"ACCEPTANCE RATIO FOR EPSILON=                  " << ar_e << "/" << num_steps << "=" <<ar_e / (double) num_steps << "\n";
}

void Arguments::print_ar_my_update(const unsigned int &num_steps){
  cout<<"ACCEPTANCE RATIO FOR LAMBDA=                    " << ar_l << "/" << K * N * num_steps  << "=" << ar_l / (double) (K * N * num_steps) << "\n";
  cout <<"ACCEPTANCE RATIO FOR EPSILON=                  " << ar_e << "/" << num_steps << "=" <<ar_e / (double) num_steps << "\n";
}

Arguments::Arguments(int argc, const char* argv[]){  
  if (argc != 31){
    cerr << argv[0] << ": incorrect number of arguments. Aborting.\n";
    exit(1);
  }
  VERBOSE = atoi(argv[1]);
  C = f_tokenize(argv[2], ",", 0);
  N = (unsigned int) C.size();
  Csum = 0.;
  for (unsigned int i = 0; i < N; i ++)
    Csum += C[i];
  NORMALIZATION = ((double) Csum) /(double) N;
  for (unsigned int i = 0; i < N; i ++){
    C.at(i) = C[i] / NORMALIZATION;
    log_C.push_back(log(C[i]));
  }
  RL = atoi(argv[3]);
  argv4 = argv[4];
  listCoordinate = f_tokenize(argv4, ",", 0);
  lC = (unsigned int) listCoordinate.size();
  MC_STEPS = atoi(argv[9]);
  MC_BURNIN = atoi(argv[10]);
  MC_EQ = atoi(argv[11]);
  MC_THIN = atoi(argv[12]);
  genePrediction = (unsigned int) atoi(argv[13]);
  genePredFile = argv[14];
  argv15 = argv[15];
  chr = argv[16];
  MIN_EX_LEN = atoi(argv[17]);
  MAX_EX_LEN_DIV_2 = (unsigned int) ceil(atoi(argv[18])/2.);//interval around a center used to choose SS
  MIN_IN_LEN = atoi(argv[19]);
  MAX_IN_LEN = atoi(argv[20]);
  scale_junction_count_by = atoi(argv[21]);

  //set parameters of the priors on log_lambda and log_epsilon                 
  mu_e = atof(argv[22]);
  var_e = atof(argv[23]);
  sd_e = pow(var_e, 0.5);
  sd_q = atof(argv[24]);
  var_l = atof(argv[25]); ///?????????????
  c = atof(argv[26]);
  a = atof(argv[27]); // var_log_lambda\sim InvGamma(a,b)
  b = atof(argv[28]);

  utils = argv[29];
  OVERHANG = atoi(argv[30]);

  LC = listCoordinate[lC - 1] - listCoordinate[0] + 1.;//length of the region       
  sumC = 0; for (unsigned int i = 0; i < N; i ++) sumC += C[i];
  //vector with the length of the intervals defined by listCoordinate    
  for (int i = 0; i < lC - 1; i ++) L.push_back((unsigned int) (listCoordinate[i + 1] - listCoordinate[i]));
  alphaSS = 1.;//dirichlet prior parameter for centers
  if (argv[5] != "NA" && argv[6] != "NA"){
    pos3 = i_tokenize(argv[5], ",", 1);
    pos5 = i_tokenize(argv[6], ",", 1);                               
  }
  if (argv[7] != "NA" && argv[8] != "NA"){
    neg3 = i_tokenize(argv[7], ",", 1);
    neg5 = i_tokenize(argv[8], ",", 1);
  }
  if (genePrediction == 1 || genePrediction == 0){
    genePredOut = argv15;
    genePredOut += "GenePredOut";
  }else 
    if (genePrediction == 2) 
      genePredOut = genePredFile;
  if (VERBOSE==1){
    outfilellk = argv15;
    outfilellk += "llkOut";
    outdatallk.open(outfilellk.c_str());
  }
  summary = argv15 + "summary.txt";
  ExprOut = argv15 + "ExprOut";
  //Initialize random seed                                             
  rSeed = random_seed();
  gsl_rng_env_setup();
  Tor = gsl_rng_default;
  rnd = gsl_rng_alloc(Tor);
  gsl_rng_set(rnd, rSeed);
  my_types.push_back("Propose Change 5SS"); 
  my_types.push_back("Propose Change 3SS"); 
  my_types.push_back("Propose Add Center"); 
  my_types.push_back("Propose Remove Center"); 
  my_types.push_back("Propose Swap"); 
  my_types.push_back("Propose Recombine");
}

void Arguments::set0(unsigned int &my_pK, unsigned int &my_nK){
  pK = my_pK; 
  nK = my_nK;
  K = pK + nK;
  powK = (unsigned int) pow(2., (double) K);
  Z = gsl_matrix_alloc(powK, K);
  for (unsigned long z = 0; z < powK; z ++){
    boost::dynamic_bitset<> ZZ(K, z);
    for (unsigned int k = 0; k < K; k ++)
      gsl_matrix_set(Z, z, k, ZZ[K - k - 1]);
  }
  //////////////////////////////////////                            
  VAR_q = gsl_matrix_alloc(K, K);//variance of the proposal on log_lambda         
  VAR_l = gsl_matrix_alloc(K, K);//variance of the prior on log_lambda       
  gsl_matrix_set_identity(VAR_q);
  gsl_matrix_scale(VAR_q, sd_q * sd_q);
  gsl_matrix_set_identity(VAR_l);
  gsl_matrix_scale(VAR_l, var_l);
}
 
void Arguments::set2(unsigned int &my_pK, unsigned int &my_nK){
  gsl_matrix_free(Z);
  gsl_matrix_free(VAR_q);
  gsl_matrix_free(VAR_l);

  set0(my_pK, my_nK);
}

