/** \file utils.h
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
