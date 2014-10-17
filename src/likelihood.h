/** \file likelihood.h
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

#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include "Arguments.h"
#include "Reads.h"

double diff_lln(const gsl_vector* log_lambda, const gsl_vector* log_gamma, const double &variance, unsigned int &K);
double diff_lln(const gsl_vector* log_lambda, const gsl_vector* log_gamma, const gsl_vector* log_lambda_bar, const gsl_vector* Var_l, Arguments &A);
double diff_lln(const gsl_vector* log_lambda, const gsl_vector* log_gamma, const gsl_vector* log_lambda_bar, const gsl_vector* Var_l, unsigned int &k, unsigned int &i, Arguments &A);
double log_likelihood(const gsl_vector *lambda, const double &epsilon, const double &log_epsilon, const gsl_matrix* X, const gsl_vector* L, allReads &reads, Arguments &A);

#endif
