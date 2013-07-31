#ifndef EXPRESSION_H
#define EXPRESSION_H

#include "TranscriptModel.h"
#include "likelihood.h"

void init_log_lambda(gsl_vector* log_lambda, Arguments &A);
void init_log_lambda_v(gsl_vector* log_lambda, Arguments &A);
void update_lambda(gsl_vector* lambda, gsl_vector* log_lambda, double &epsilon, double &log_epsilon, TranscriptModel &transcriptModel, double &ll, Arguments &A);
void update_epsilon (gsl_vector* lambda, double &epsilon, double &log_epsilon, TranscriptModel &transcriptModel, double &ll, Arguments &A);

void print_to_Abuffer(gsl_vector* lambda, const double &epsilon, const double &ll, Arguments &A);

#endif
