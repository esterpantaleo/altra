#include "expression.h"

void print_to_Abuffer(gsl_vector* lambda, const double &epsilon, const double &ll, Arguments &A){
  for (unsigned int i = 0; i < A.N; i ++){
    for (unsigned int k = 0; k < A.K; k ++){
      A.buffer << setprecision(4) << gsl_vector_get(lambda, A.K * i + k) << " ";
    }
  }
  A.buffer << setprecision(4) << epsilon << " " << fixed << ll << "\n";
}

//N=1                                                                          
void init_log_lambda(gsl_vector* log_lambda, Arguments &A){
  //initialize log_lambda                                           
  rmvnorm(A.rnd, A.VAR_l, log_lambda);
}

void init_log_lambda_v(gsl_vector* log_lambda, Arguments &A){
  double V_l, sqrtc = pow(A.c, 0.5);
  gsl_vector* log_lambda_bar = gsl_vector_alloc(A.K);
  gsl_vector* Var_l = gsl_vector_alloc(A.K);

  for (unsigned int k = 0; k < A.K; k ++){
    V_l = 1. / gsl_ran_gamma(A.rnd, A.a, 1. / A.b);
    gsl_vector_set(Var_l, k, V_l);
    V_l = pow(V_l, 0.5);
    gsl_vector_set(log_lambda_bar, k, gsl_ran_gaussian(A.rnd, V_l * sqrtc));
    for (unsigned int i = 0; i < A.N; i ++)//scalar product...                  
      gsl_vector_set(log_lambda, A.K * i + k, gsl_vector_get(log_lambda_bar, k) + gsl_ran_gaussian(A.rnd, V_l));
  }

  gsl_vector_free(log_lambda_bar);
  gsl_vector_free(Var_l);
}

void update_lambda(gsl_vector* lambda, gsl_vector* log_lambda, double &epsilon, double &log_epsilon, TranscriptModel &transcriptModel, double &ll, Arguments &A){
  double f = 1. + (double) A.N * A.c, d = A.c / f, alpha, m, s2, v, tmp, sn2, m2, llproposed;
  gsl_vector* log_gamma = gsl_vector_alloc(A.N * A.K);
  gsl_vector* gamma = gsl_vector_alloc(A.N * A.K);
  gsl_vector* log_lambda_bar = gsl_vector_alloc(A.K);
  gsl_vector* Var_l = gsl_vector_alloc(A.K);

  //draw Var_l and log_lambda_bar 
  for (unsigned int k = 0; k < A.K; k ++){
    m = 0.;
    s2 = 0;
    for (unsigned int i = 0; i < A.N; i ++)
      m += gsl_vector_get(log_lambda, A.K * i + k);
    m /= A.N;
    for (unsigned int i = 0; i < A.N; i ++)
      s2 += (gsl_vector_get(log_lambda, A.K * i + k) - m) * (gsl_vector_get(log_lambda, A.K * i + k) - m);
    m2 = m * m;
    sn2 = A.b + 0.5 * (s2 + m2 * A.N / f);
    v = 1. / gsl_ran_gamma(A.rnd, A.a + 0.5 * ((double) A.N), 1. / sn2);
    //update variance log_lambda and log_lambda_bar                                                      
    gsl_vector_set(Var_l, k, v);
    gsl_vector_set(log_lambda_bar, k, d * m * A.N + gsl_ran_gaussian(A.rnd, pow(d * v, 0.5)));
  }
  //propose log_gamma    
  //compute alpha                                                                                               
  for (unsigned int k = 0; k < A.K; k ++){
    for (unsigned int i = 0; i < A.N; i ++){
      gsl_vector_memcpy(log_gamma, log_lambda);
      gsl_vector_set(log_gamma, A.K * i + k, gsl_vector_get(log_gamma, A.K * i + k) + gsl_ran_gaussian(A.rnd, A.sd_q));
      gsl_vector_memcpy(gamma, lambda);
      gsl_vector_set(gamma, A.K * i + k, exp(gsl_vector_get(log_gamma, A.K * i + k)));

      llproposed = log_likelihood(gamma, epsilon, log_epsilon, transcriptModel.X, transcriptModel.L, A);
      alpha = llproposed - ll;
      alpha += diff_lln(log_lambda, log_gamma, log_lambda_bar, Var_l, k, i, A);
      alpha = min(0., alpha);
      if (gsl_rng_uniform(A.rnd) < exp(alpha)){
	A.ar_l += 1;//accept                                      
	gsl_vector_swap(log_lambda, log_gamma);
	gsl_vector_swap(lambda, gamma);
	ll = llproposed;
      }
    }
  }

  //Free                                                            
  gsl_vector_free(log_gamma);
  gsl_vector_free(gamma);
  gsl_vector_free(log_lambda_bar);
  gsl_vector_free(Var_l);

  return;
}

void update_epsilon (gsl_vector* lambda, double &epsilon, double &log_epsilon, TranscriptModel &transcriptModel, double &ll, Arguments &A){
  double log_delta, delta, alpha, llproposed;
  log_delta = log_epsilon + gsl_ran_gaussian(A.rnd, A.sd_q);
  delta = exp(log_delta);
  //calculate alpha                                           
  llproposed = log_likelihood(lambda, delta, log_delta, transcriptModel.X, transcriptModel.L, A);
  alpha = llproposed - ll;
  alpha += 0.5 * (log_epsilon * log_epsilon - log_delta * log_delta - 2. * A.mu_e * (log_epsilon - log_delta)) / A.var_e;
  alpha = min(0., alpha);
  if (gsl_rng_uniform(A.rnd) < exp(alpha)){
    A.ar_e ++;//accept                                  
    log_epsilon = log_delta;
    epsilon = delta;
    ll = llproposed;
  }//else reject and stay in log_epsilon        
}
