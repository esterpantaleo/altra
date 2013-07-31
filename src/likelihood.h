#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include "Arguments.h"

double diff_lln(const gsl_vector* log_lambda, const gsl_vector* log_gamma, const double &variance, unsigned int &K);
double diff_lln(const gsl_vector* log_lambda, const gsl_vector* log_gamma, const gsl_vector* log_lambda_bar, const gsl_vector* Var_l, Arguments &A);
double diff_lln(const gsl_vector* log_lambda, const gsl_vector* log_gamma, const gsl_vector* log_lambda_bar, const gsl_vector* Var_l, unsigned int &k, unsigned int &i, Arguments &A);
double log_likelihood(const gsl_vector *lambda, const double &epsilon, const double &log_epsilon, const gsl_matrix* X, const gsl_vector* L, Arguments &A);

#endif
