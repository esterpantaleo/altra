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
  ~Arguments(){
    outdatallk.close();
    outdata1.close();
    gsl_matrix_free(Z);
    gsl_matrix_free(VAR_q);
    gsl_matrix_free(VAR_l);
    gsl_rng_free(rnd);
  }
  void set0(unsigned int &my_pK, unsigned int &my_nK);
  void set(unsigned int &my_pK, unsigned int &my_nK);
  void set2(unsigned int &my_pK, unsigned int &my_nK);
  
  void clear_ar(){
    ar_l = 0.;
    ar_e = 0.;
  }

  void print_ar(const unsigned int &num_steps){
    cout<<"ACCEPTANCE RATIO FOR LAMBDA=                    " << ar_l << "/" << N * num_steps  << "=" << ar_l / (double) (N * num_steps) << "\n";
    cout <<"ACCEPTANCE RATIO FOR EPSILON=                  " << ar_e << "/" << num_steps << "=" <<ar_e / (double) num_steps << "\n";
  }
  void print_ar_my_update(const unsigned int &num_steps){
    cout<<"ACCEPTANCE RATIO FOR LAMBDA=                    " << ar_l << "/" << K * N * num_steps  << "=" << ar_l / (double) (K * N * num_steps) << "\n";
    cout <<"ACCEPTANCE RATIO FOR EPSILON=                  " << ar_e << "/" << num_steps << "=" <<ar_e / (double) num_steps << "\n";
  }
};

Arguments::Arguments(int argc, const char* argv[]){  
  if (argc != 30){
    cerr << argv[0] << ": incorrect number of arguments. Aborting.\n";
    exit(1);
  }
  C = f_tokenize(argv[1], ",", 0);
  N = (unsigned int) C.size();
  Csum = 0.;
  for (unsigned int i = 0; i < N; i ++)
    Csum += C[i];
  NORMALIZATION = ((double) Csum) /(double) N;
  for (unsigned int i = 0; i < N; i ++){
    C.at(i) = C[i] / NORMALIZATION;
    log_C.push_back(log(C[i]));
  }
  RL = atoi(argv[2]);
  argv3 = argv[3];
  listCoordinate = f_tokenize(argv3, ",", 0);
  lC = (unsigned int) listCoordinate.size();
  MC_STEPS = atoi(argv[8]);
  MC_BURNIN = atoi(argv[9]);
  MC_EQ = atoi(argv[10]);
  MC_THIN = atoi(argv[11]);
  genePrediction = (unsigned int) atoi(argv[12]);
  genePredFile = argv[13];
  argv14 = argv[14];
  chr = argv[15];
  MIN_EX_LEN = atoi(argv[16]);
  MAX_EX_LEN_DIV_2 = (unsigned int) ceil(atoi(argv[17])/2.);//interval around a center used to choose SS
  MIN_IN_LEN = atoi(argv[18]);
  MAX_IN_LEN = atoi(argv[19]);
  scale_junction_count_by = atoi(argv[20]);

  //set parameters of the priors on log_lambda and log_epsilon                 
  mu_e = atof(argv[21]);
  var_e = atof(argv[22]);
  sd_e = pow(var_e, 0.5);
  sd_q = atof(argv[23]);
  var_l = atof(argv[24]); ///?????????????
  c = atof(argv[25]);
  a = atof(argv[26]); // var_log_lambda\sim InvGamma(a,b)
  b = atof(argv[27]);

  utils = argv[28];
  output2summary = argv[29];

  LC = listCoordinate[lC - 1] - listCoordinate[0] + 1.;//length of the region       
  sumC = 0; for (unsigned int i = 0; i < N; i ++) sumC += C[i];
  //vector with the length of the intervals defined by listCoordinate    
  for (int i = 0; i < lC - 1; i ++) L.push_back((unsigned int) (listCoordinate[i + 1] - listCoordinate[i]));
  alphaSS = 1.;//dirichlet prior parameter for centers
  if (argv[4] != "NA" && argv[5] != "NA"){
    pos3 = i_tokenize(argv[4], ",", 1);
    pos5 = i_tokenize(argv[5], ",", 1);                               
  }
  if (argv[6] != "NA" && argv[7] != "NA"){
    neg3 = i_tokenize(argv[6], ",", 1);
    neg5 = i_tokenize(argv[7], ",", 1);
  }
  if (genePrediction == 1 || genePrediction == 0){
    genePredOut = argv14;
    genePredOut += "GenePredOut";
  }else 
    if (genePrediction == 2) 
      genePredOut = genePredFile;
  
  outfilellk = argv14;
  outfilellk += "llkOut";
  
  outdatallk.open(outfilellk.c_str());
  outfile1 = argv14;
  outfile1 += "MC_output1";
  outdata1.open(outfile1.c_str());
  
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

 
void Arguments::set(unsigned int &my_pK, unsigned int &my_nK){
  set0(my_pK, my_nK);
  unsigned int PRINTED_STEPS = (unsigned int) ((double)MC_STEPS/(double)MC_THIN) - 1;
  
  run_output2summary = "source \"";
  run_output2summary += utils;
  run_output2summary += "\";";
  run_output2summary += output2summary;
  run_output2summary += " ";
  run_output2summary += argv14;
  run_output2summary += " ";
  run_output2summary += ItoA(PRINTED_STEPS);
  run_output2summary += " ";
  run_output2summary += ItoA(pK);
  run_output2summary += " ";
  run_output2summary += ItoA(nK);
  run_output2summary += " ";
  run_output2summary += chr;
  run_output2summary += " ";
  run_output2summary += argv3;
}

void Arguments::set2(unsigned int &my_pK, unsigned int &my_nK){
  gsl_matrix_free(Z);
  gsl_matrix_free(VAR_q);
  gsl_matrix_free(VAR_l);

  set0(my_pK, my_nK);
}
    
#endif //ARGUMENTS_H 
