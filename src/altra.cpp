/** \file altra.cpp
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


#include "expression.h"

vector<TypeOfMove> Move::move_proposed = vector<TypeOfMove>();
vector<ModifiedCenter> Move::modifiedCenter = vector<ModifiedCenter>();
vector<unsigned int> Move::t = vector<unsigned int>();
vector<SS*> Move::s = vector<SS*>(); 
Centers* Move::centers = NULL;

int Move::counter = 0;
int Move::ar_counter = 0;
int Move::small_world = 0;

int main (int argc, const char* argv[]){
  unsigned int step;
  
  Arguments A(argc, argv);
  allReads reads(A);
  double log_epsilon, epsilon, ll;
  gsl_vector* lambda;
  gsl_vector* log_lambda;

  if (A.genePrediction == 1){
    //**************************************************************************************************     
    //    INITIALIZE TRANSCRIPT MODEL, i.e., posCenters, negCenters, X and L           
    //**************************************************************************************************     
    TranscriptModel transcriptModel(A.genePredFile, reads, A);

    //**************************************************************************************************
    //     INITIALIZE LAMBDAs AND EPSILON                          
    //**************************************************************************************************  
    //initialize lambda and epsilon
    lambda = gsl_vector_alloc(A.N * A.K);
    log_lambda = gsl_vector_alloc(A.N * A.K);
    if (A.N == 1){
      for (unsigned int i = 0; i < A.N; i ++)
	init_log_lambda(&gsl_vector_subvector(log_lambda, A.K * i, A.K).vector, A);
    }else
      init_log_lambda_v(log_lambda, A);
    exp(lambda, log_lambda, A.K * A.N);
    log_epsilon = A.mu_e + gsl_ran_gaussian(A.rnd, A.sd_e);
    epsilon = exp(log_epsilon);
  
    //print initial value
    cout<<"initial state:\n"; transcriptModel.print_to_Abuffer(ll,A);  cout << A.buffer.str();  A.buffer.str("");
    cout << "X="; transcriptModel.print_X(); cout << "L = "; transcriptModel.print_L();
    ll = log_likelihood(lambda, epsilon, log_epsilon, transcriptModel.X, transcriptModel.L, A);
    print_to_Abuffer(lambda, epsilon, ll, A); cout<<A.buffer.str(); A.buffer.str(""); 

    //*****************************************************************************************************
    //     RUN MC; print to outdata
    //*****************************************************************************************************    
    string outfile = A.argv14;
    outfile += "MC_output";
    ofstream outdata;
    outdata.open(outfile.c_str());
    A.clear_ar();
    
    for (step = 0; step < A.MC_STEPS; step ++){
      //if (step % 1000 == 0 && step != 0) cout<<"step "<<step<<"...\n";     
      //update transcript
      transcriptModel.update(reads, lambda, epsilon, log_epsilon, ll, A);
      update_lambda(lambda, log_lambda, epsilon, log_epsilon, transcriptModel, ll, A);
      update_epsilon(lambda, epsilon, log_epsilon, transcriptModel, ll, A);
    
      if (step % A.MC_THIN == 0){
	transcriptModel.print_to_Abuffer(ll,A); outdata << A.buffer.str(); A.buffer.str("");
      }
      if (step % 10 == 0) A.outdatallk << ll << "\n";
    }
    
    //print final state 
    A.print_ar_my_update(step); transcriptModel.print_ar(A);
    cout<<"ar_small_world="<<Move::ar_counter<<"/"<<Move::counter<<"="<< (double) Move::ar_counter / (double) Move::counter <<"\n";
    cout<<"last final state:\n";cout << "X="; transcriptModel.print_X(); cout << "L = "; transcriptModel.print_L();
    transcriptModel.print_to_Abuffer(ll,A); cout << A.buffer.str();A.buffer.str("");
    print_to_Abuffer(lambda, epsilon, ll, A); cout<<A.buffer.str();A.buffer.str("");
    //free                          
    outdata.close();
    
    //****************************************************************************************************
    //   SUMMARIZE_RESULTS
    //****************************************************************************************************
    cout << A.run_output2summary.c_str()<<"\n"; system(A.run_output2summary.c_str());

    gsl_vector_free(lambda);
    gsl_vector_free(log_lambda);
  }

  //scale back to reads with junctions                                          
  if (A.scale_junction_count_by != 1){
    for (unsigned int i = 0; i < A.N; i ++){
      for (unsigned int j = 0; j < reads.pos.size; j ++)
	reads.pos.v[j].counts[i] /= A.scale_junction_count_by;
      for (unsigned int j = 0; j < reads.neg.size; j ++)
	reads.neg.v[j].counts[i] /= A.scale_junction_count_by;
    }
  }
  
  TranscriptModel optimalTranscriptModel(A.genePredOut, reads, A);//define the new A.K A.pK A.nK 

  lambda = gsl_vector_alloc(A.N * A.K);
  log_lambda = gsl_vector_alloc(A.N * A.K);
  //initialize lambda and epsilon                                
  if (A.N == 1)
    init_log_lambda(log_lambda, A);
  else
    init_log_lambda_v(log_lambda, A);
  exp(lambda, log_lambda, A.K * A.N);
  log_epsilon = A.mu_e + gsl_ran_gaussian(A.rnd, A.sd_e);
  epsilon = exp(log_epsilon);
  ll = log_likelihood(lambda, epsilon, log_epsilon, optimalTranscriptModel.X, optimalTranscriptModel.L, A);

  //**************************************************************************************************     
  //    FIND OPTIMAL LAMBDAS AND EPSILON FOR FINAL STATE; print to outdata1                  
  //************************************************************************************************** 
  A.clear_ar();
 
  for (step = 0; step < A.MC_EQ + A.MC_BURNIN; step ++){
    update_lambda(lambda, log_lambda, epsilon, log_epsilon, optimalTranscriptModel, ll, A);
    update_epsilon(lambda, epsilon, log_epsilon, optimalTranscriptModel, ll, A);
   
    //print
    if (step > A.MC_BURNIN){
      if (step % A.MC_THIN == 0){
	print_to_Abuffer(lambda, epsilon, ll, A); A.outdata1 << A.buffer.str(); A.buffer.str("");
      }
    }
    if (step % 10 == 0) A.outdatallk << ll << "\n";
  }

  //print final state
  A.print_ar_my_update(step); cout<<"final state:"; optimalTranscriptModel.print_to_Abuffer(ll, A); cout<<"buffer="<<A.buffer.str()<<"\n";    
  print_to_Abuffer(lambda, epsilon, ll, A); cout<<A.buffer.str(); A.buffer.str("");
  cout << "X="; optimalTranscriptModel.print_X(); cout << "L = "; optimalTranscriptModel.print_L();
   
  //free
  gsl_vector_free(lambda);
  gsl_vector_free(log_lambda);
  
  return 0;
}
