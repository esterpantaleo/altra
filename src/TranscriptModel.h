#ifndef TRANSCRIPTMODEL_H
#define TRANSCRIPTMODEL_H

#include "ProposedTranscript.h"
#include "Reads.h"
#include "Move.h"
#include "likelihood.h"

//check that genePrediction is either 0 or 1 or 2
class TranscriptModel{
public:
  unsigned int K;
  Centers posCenters, negCenters;
  gsl_matrix *X;
  gsl_vector *L;
  gsl_matrix *Xproposed;
  gsl_vector *Lproposed;
  vector<boost::dynamic_bitset<> > pZ, nZ, uZ, pZproposed, nZproposed, uZproposed, pZproposed2, nZproposed2, uZproposed2;
  vector<ProposedTranscript> proposedTranscripts;
   
  void print_ar(Arguments &A);
  void setZ (allReads &reads, Arguments &A);
  void print_Z(vector<boost::dynamic_bitset<> > &aZ);
  void setX(gsl_matrix* my_X, Reads &reads, vector<boost::dynamic_bitset<> > &sZ, Arguments &A);
  void setX(allReads &reads, Arguments &A);
  void setL(Arguments &A);
  void setZproposed(Reads &uns, Reads &my_pos, vector<boost::dynamic_bitset<> > &my_pZ, vector<boost::dynamic_bitset<> > &my_pZproposed, unsigned int &AK, Arguments &A);
  void setXproposed(allReads &reads, Reads &my_pos, Reads &my_neg, vector<boost::dynamic_bitset<> > &my_pZ, vector<boost::dynamic_bitset<> > &my_pZproposed, vector<boost::dynamic_bitset<> > &my_nZ, unsigned int AK, Arguments &A);
  void setLproposed(unsigned int AK, Arguments &A);

  TranscriptModel(string &genePredFile, allReads &reads, Arguments &A);
  TranscriptModel(string &genePredFile, Centers &old_posCenters, Centers &old_negCenters, allReads &reads, Arguments &A);
  ~TranscriptModel();

  static Move *draw_move(ModifiedCenter &modifiedCenter, Arguments &A);
  void reject();
  void accept_reject(allReads &reads, Reads &my_pos, Reads &my_neg, vector<boost::dynamic_bitset<> > &my_pZ, vector<boost::dynamic_bitset<> > &my_pZproposed, vector<boost::dynamic_bitset<> > &my_nZ, unsigned int AK, gsl_vector* lambda, double &epsilon, double &log_epsilon, double &ll, Arguments &A); 
  void accept_reject_pos(allReads &reads, gsl_vector* lambda, double &epsilon, double &log_epsilon, double &ll, Arguments &A);
  void accept_reject_neg(allReads &reads, gsl_vector *lambda, double &epsilon, double &log_epsilon, double &ll, Arguments &A);
  void propose_update(unsigned int &t1, allReads &reads, gsl_vector *lambda, double &epsilon, double &log_epsilon, double &ll, unsigned int &K, Arguments &A);
  void update(allReads &reads, gsl_vector *lambda, double &epsilon, double &log_epsilon, double &ll, Arguments &A);
  void print_X();
  void print_L();
  void print(string strand, Centers &centers, Arguments &A);
  void print_to_Abuffer(const double &ll, Arguments &A);
};

#endif
