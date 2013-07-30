#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <math.h>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <functional>
//#include <stdlib.h>
//#include <gzstream.h>
#include <utility>
#include <gsl/gsl_rng.h>
//#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_multimin.h>
#include <boost/dynamic_bitset.hpp>
#include <iomanip>

extern "C"
{
#include <assert.h>
  // LU decomoposition of a general matrix                                                                                                                                                                  
  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  // generate inverse of a matrix given its LU decomposition                                                                                                                                                
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}


using namespace std;

vector<string> string_tokenize(const string& str, const string& delimiters, bool skip_empty);
vector<unsigned int> i_tokenize(const string& str, const string& delimiters, bool skip_empty);
vector<double> f_tokenize(const string& str, const string& delimiters, bool skip_empty);
unsigned long int random_seed();
void rmvnorm(const gsl_rng *r, const gsl_matrix *V, gsl_vector *r_sample);
string ItoA(unsigned int &my_int);
void exp(gsl_vector *gamma, const gsl_vector *log_gamma, unsigned int size);
void print_my_gsl_vector(gsl_vector *v);
void print_my_gsl_matrix(gsl_matrix *m);

#endif //UTILS_H 
