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
#include "Summary.h"

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
  
    ll = log_likelihood(lambda, epsilon, log_epsilon, transcriptModel.X, transcriptModel.L, reads, A);

    //print initial state
    if (A.VERBOSE==1){
       // print loglikelihood
       cout<<"initial state:\n"; 
       transcriptModel.print_to_Abuffer(ll,A);  
       cout << A.buffer.str();  
       A.buffer.str("");
       // print X and L
       cout << "X="; 
       transcriptModel.print_X(); 
       cout << "L="; 
       transcriptModel.print_L();
       // print lambda and epsilon
       print_to_Abuffer(lambda, epsilon, ll, A); 
       cout<<A.buffer.str(); 
       A.buffer.str(""); 
    }

    //*****************************************************************************************************
    //     RUN MC
    //*****************************************************************************************************    
    A.clear_ar();

    //initialise a map
    //Create a c++ map of the explored states
    //The _key_ of the map is transcriptModel.z (the bitset);
    //the _value_ of the map is _a_vector_of_likelihoods_
    //containing all the likelihoods of the system when it was exploring _key_
    //_length(value)_ is then the number of times the system explored _key_ 
    SummaryAsMap myMap;  
    for (step = 0; step < A.MC_STEPS; step ++){
      //if (step % 1000 == 0 && step != 0) cout<<"step "<<step<<"...\n";     
      //update transcript model
      transcriptModel.update(reads, lambda, epsilon, log_epsilon, ll, A);
      update_lambda(lambda, log_lambda, epsilon, log_epsilon, transcriptModel, ll, reads, A);
      update_epsilon(lambda, epsilon, log_epsilon, transcriptModel, ll, reads, A);

      //push_back to map
      if (step % A.MC_THIN == 0 && step > A.MC_BURNIN){
        //convert the TranscriptModel (e.g.:100110,110100,101010)
        //into a bitset (e.g.:100110110100101010)
        transcriptModel.to_z(A);
        myMap[transcriptModel.z].push_back(ll);
      }
      if (A.VERBOSE==1){
	//print loglikelihood
         if (step % A.MC_THIN == 0) A.outdatallk << ll << "\n";
      }
    }
    
    //print final state
    if (A.VERBOSE==1){
      //print ar (acceptance ratios)
       A.print_ar_my_update(step); 
       transcriptModel.print_ar(A);
       cout<<"ar_small_world="<<Move::ar_counter<<"/"<<Move::counter<<"="<< (double) Move::ar_counter / (double) Move::counter <<"\n";
       //print final state X and L
       cout<<"last state:\n";
       cout << "X="; transcriptModel.print_X(); 
       cout << "L="; transcriptModel.print_L();
       //print transcript model and loglikelihood
       transcriptModel.print_to_Abuffer(ll,A); 
       cout << A.buffer.str();
       A.buffer.str("");
       //print lambda and epsilon
       print_to_Abuffer(lambda, epsilon, ll, A); 
       cout<<A.buffer.str();
       A.buffer.str("");
    }

    //****************************************************************************************************
    //   SUMMARIZE_RESULTS
    //****************************************************************************************************
    //summarize myMap elements into a Summary that can then be sorted
    //based on the average loglkelihood  
    Summary summary(&myMap);
    myMap.clear();  
    summary.print(A); //print summary and genepred of first line
  
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
  ll = log_likelihood(lambda, epsilon, log_epsilon, optimalTranscriptModel.X, optimalTranscriptModel.L, reads, A);

  //**************************************************************************************************     
  //    FIND OPTIMAL LAMBDAS AND EPSILON FOR FINAL STATE
  //************************************************************************************************** 
  A.clear_ar();

  vector<vector<double> > lambdas;
  vector<double> epsilons;
  for (step = 0; step < A.MC_EQ + A.MC_BURNIN; step ++){
    update_lambda(lambda, log_lambda, epsilon, log_epsilon, optimalTranscriptModel, ll, reads, A);
    update_epsilon(lambda, epsilon, log_epsilon, optimalTranscriptModel, ll, reads, A);
   
    //MC step
    if (step > A.MC_BURNIN)
      if (step % A.MC_THIN == 0){
	vector<double> lambda_vector;
	for (unsigned int i = 0; i < A.N; i ++)   
	  for (unsigned int k = 0; k < A.K; k ++)
	    lambda_vector.push_back(gsl_vector_get(lambda, A.K * i + k));
	lambdas.push_back(lambda_vector);
	epsilons.push_back(epsilon);
	//clean
	lambda_vector.clear();
      
	//print loglikelihood
	if (A.VERBOSE)
	  A.outdatallk << ll << "\n";
      }
  }

  //print final state
  if (A.VERBOSE==1){
    //print ar
    A.print_ar_my_update(step); 
    //print transcript model and loglikelihood 
    cout << "final state:" << endl; 
    optimalTranscriptModel.print_to_Abuffer(ll, A); 
    cout << A.buffer.str();
    A.buffer.str("");
    //print lambda and epsilon
    cout << "lambdas, epsilon, llk:" <<endl;    
    print_to_Abuffer(lambda, epsilon, ll, A); 
    cout << A.buffer.str(); 
    A.buffer.str("");
    //print final state X and L
    cout << "X="; 
    optimalTranscriptModel.print_X(); 
    cout << "L="; 
    optimalTranscriptModel.print_L();
  }
  unsigned int nsteps = (unsigned int) epsilons.size();
  double epsilon_sum = 0, epsilon_sumsq = 0;
  vector<double> lambda_sum(A.N*A.K, 0.0), lambda_sumsq(A.N*A.K, 0.0);
  for (unsigned int j=0; j<nsteps; j++){
    for (unsigned int i=0; i<A.N*A.K; i++){
      lambda_sum[i] += lambdas[j][i];
      lambda_sumsq[i] += lambdas[j][i]*lambdas[j][i];
    }
  }
  for (unsigned int j=0; j<nsteps; j++){
    epsilon_sum += epsilons[j];
    epsilon_sumsq += epsilons[j]*epsilons[j];
  }
  //////////////print expression values to file 
  ofstream ssfile;
  ssfile.open(A.ExprOut.c_str());
  //print lambda
  for (unsigned int i = 0; i < A.N; i ++){
    for (unsigned int k = 0; k < A.K; k ++){
      ssfile << lambda_sum[i*A.K+k]/nsteps << " " << sqrt(lambda_sumsq[i*A.K+k]/nsteps - (lambda_sum[i*A.K+k]*lambda_sum[i*A.K+k])/(nsteps*nsteps)) << " ";
    }
  }
  //print epsilon
  ssfile << epsilon_sum/nsteps << " " << sqrt(epsilon_sumsq/nsteps - epsilon_sum*epsilon_sum/(nsteps*nsteps)) << endl;
  ssfile.close();

  //free
  gsl_vector_free(lambda);
  gsl_vector_free(log_lambda);
  
  return 0;
}
