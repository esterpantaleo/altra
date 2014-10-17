/** \file utils.cpp
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

#include "utils.h"

//tokenize a string into a vector of strings                                                    
vector<string> string_tokenize(const string& str, const string& delimiters, bool skip_empty) {
  // Skip delimiters at beginning.                                                                                                                  
  string::size_type lastPos = skip_empty ? str.find_first_not_of(delimiters, 0) : 0;
  // Find first "non-delimiter".                                                                                                                                                                           
  string::size_type pos=str.find_first_of(delimiters, lastPos);
  vector<string> toreturn;
  toreturn.clear();

  while (string::npos != pos || string::npos != lastPos){
    // Found a token, add it to the vector.                                                                                                                                                                
    //__ASSERT(pos > lastPos || !skip_empty, "internal error, pos <= lastPos.\n");                                                                                                                         
    //if (pos == lastPos) toreturn.push_back("");                                                                                                                                                          
    toreturn.push_back(str.substr(lastPos, pos - lastPos));
    if (pos == string::npos) break;
    if (pos == str.length() - 1) {
      if (!skip_empty) toreturn.push_back("");
      break;
    }
    // Skip delimiters.  Note the "not_of"                                                            
    lastPos = skip_empty ? str.find_first_not_of(delimiters, pos) : pos + 1;
    // Find next "non-delimiter"                                               
    pos = str.find_first_of(delimiters, lastPos);
  }

  return toreturn;
}

//tokenize a string into a vector of int                                                
vector<unsigned int> i_tokenize(const string& str, const string& delimiters, bool skip_empty) {
  vector<string> tokens=string_tokenize(str, delimiters, skip_empty);
  //int I=(int) tokens.size();                              
  vector<unsigned int> toreturn;
  //if (I>0){                                                                                                                            
  for (int i=0;i<tokens.size();i++)
    toreturn.push_back((unsigned int)atoi(tokens[i].c_str()));
  //}                                                                                                                 

  return toreturn;
}

//tokenize a string into a vector of double                                                          
vector<double> f_tokenize(const string& str, const string& delimiters, bool skip_empty) {
  vector<string> tokens=string_tokenize(str, delimiters, skip_empty);
  //int I=(int) tokens.size();                                                                                      
  vector<double> toreturn;
  //if (I>0){                                                                                                                                                                                              
  for (int i=0;i<tokens.size();i++)
    toreturn.push_back(atof(tokens[i].c_str()));
  //}                                                                                                                                                 
  return toreturn;
}

unsigned long int random_seed(){
  unsigned int seed;
  FILE *devrandom;

  if ((devrandom = fopen("/dev/random","r")) == NULL) {
    fprintf(stderr,"Cannot open /dev/random, setting seed to 0\n");
    seed = 0;
  }
  else{
    fread(&seed,sizeof(seed),1,devrandom);
    fclose(devrandom);
  }

  return(seed);
}

//the following function generates a random sample from a multivariate normal distribution MVN(0,V)     
//to a vector r_sample, with V a d*d covariance matrix.                                                 
void rmvnorm(const gsl_rng *r, const gsl_matrix *V, gsl_vector *r_sample){
  //if (V->size1!=V->size2) {cout<<"the normal covariance matrix V is not a square matrix"<<endl; return;} 
  size_t d=V->size1;
  gsl_matrix *V_temp = gsl_matrix_alloc(d,d);
  gsl_matrix_memcpy(V_temp,V);
  gsl_linalg_cholesky_decomp(V_temp);
  for(size_t i=0; i<d; ++i)
    gsl_vector_set(r_sample, i, gsl_ran_ugaussian(r) );
  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, V_temp, r_sample);
  gsl_matrix_free(V_temp);

  return;
}

string ItoA(unsigned int &my_int){
  stringstream ss;
  ss << my_int;
  return ss.str();
}

//compute gamma from log_gamma                                                     
void exp(gsl_vector *gamma, const gsl_vector *log_gamma, unsigned int size){
  for (unsigned int i = 0; i < size; i ++)
    gsl_vector_set(gamma, i, exp(gsl_vector_get(log_gamma, i)));
}

void print_my_gsl_vector(gsl_vector *v){
  for (int hh = 0;hh < v -> size; hh ++)
    cout << gsl_vector_get(v, hh) << ",";
  cout << "\n";
}

void print_my_gsl_matrix(gsl_matrix *m){
  for (int hh = 0; hh < m -> size1; hh ++)
    for (int kk = 0; kk < m -> size2; kk ++)
      cout << gsl_matrix_get(m, hh, kk) << ",";
  cout << "\n";
}
