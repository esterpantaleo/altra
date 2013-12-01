/** \file likelihood.cpp
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

#include "likelihood.h"

double diff_lln(const gsl_vector* log_lambda, const gsl_vector* log_gamma, const double &variance, unsigned int &K){
  //log likelihood normal with 0 mean 
  double result, toreturn;

  gsl_blas_ddot(log_gamma, log_gamma, &result);
  toreturn = - result;
  gsl_blas_ddot(log_lambda, log_lambda, &result);
  toreturn += result;

  return toreturn /= 2. * variance;
}

//-lambda^2 + gamma^2
double diff_lln(const gsl_vector* log_lambda, const gsl_vector* log_gamma, const gsl_vector* log_lambda_bar, const gsl_vector* Var_l, Arguments &A){
  //prior ratio
  double sumsq, s, tmp, tmp2, toreturn = 0.;
  for (unsigned int k = 0; k < A.K; k ++){
    sumsq = 0.;
    s = 0.;
    for (unsigned int i = 0; i < A.N; i ++){
      //*********************************
      //to speed up use stride and ddot
      //gsl_vector_const_subvector_with_stride (const gsl_vector * v, size_t offset, size_t stride, size_t n)
      tmp = gsl_vector_get(log_lambda, A.K * i + k);
      tmp2 = gsl_vector_get(log_gamma, A.K * i + k);
      sumsq += tmp * tmp - tmp2 * tmp2;
      s += tmp2 - tmp;
    }
    toreturn += (sumsq + 2. * gsl_vector_get(log_lambda_bar, k) * s) / gsl_vector_get(Var_l, k);
  }
  return 0.5 * toreturn;
}

double diff_lln(const gsl_vector* log_lambda, const gsl_vector* log_gamma, const gsl_vector* log_lambda_bar, const gsl_vector* Var_l, unsigned int &k, unsigned int &i, Arguments &A){
  //prior ratio         
  double tmp = gsl_vector_get(log_lambda, A.K * i + k);
  double tmp2 = gsl_vector_get(log_gamma, A.K * i + k);
  return  0.5 * (tmp * tmp - tmp2 * tmp2 + 2. * gsl_vector_get(log_lambda_bar, k) * (tmp2 - tmp)) / gsl_vector_get(Var_l, k);
}

double log_likelihood(const gsl_vector *lambda, const double &epsilon, const double &log_epsilon, const gsl_matrix* X, const gsl_vector* L, Arguments &A){
  double inn[A.N], toreturn = -A.LC * epsilon * A.sumC, tmp;

  for (unsigned int i = 0; i < A.N; i ++){
    tmp = 0.;
    for (unsigned int k = 0; k < A.K; k ++)
      tmp += gsl_vector_get(L, k) * gsl_vector_get(lambda, A.K * i + k);
    toreturn -= tmp * A.C[i];
    toreturn += gsl_matrix_get(X, i, 0) * (log_epsilon + A.log_C[i]);
  }
  for (unsigned int z = 1; z < A.powK; z ++){
    for (unsigned int i = 0; i < A.N; i ++)
      inn[i] = epsilon;
    for (unsigned int k = 0; k < A.K; k ++){
      if (gsl_matrix_get(A.Z, z, k) == 1){
        for (unsigned int i = 0; i < A.N; i ++)
          inn[i] += gsl_vector_get(lambda, A.K * i + k);//modify scalar product        
      }
    }
    for (unsigned int i = 0; i < A.N; i ++)
      toreturn += gsl_matrix_get(X, i, z) * (log(inn[i]) + A.log_C[i]);
  }

  return toreturn;
}
