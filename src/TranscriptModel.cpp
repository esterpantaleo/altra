#include "TranscriptModel.h"

void TranscriptModel::print_ar(Arguments &A){
  posCenters.print_ar("+", A);
  negCenters.print_ar("-", A);
}

void TranscriptModel::setZ (allReads &reads, Arguments &A){
  uZ.clear();
  pZ.clear();
  nZ.clear();
  
  //(if center is positive only work with uns and pos; if center is negative only work with uns and neg
  for (unsigned int r = 0; r < reads.uns.size; r ++){
    boost::dynamic_bitset<> Z_r(A.K, 0);
    for (unsigned int k = 0; k < A.pK; k ++)
      Z_r[A.K - k - 1] = reads.uns.v[r].is_compatible_w(posCenters.transcripts[k], A);
    for (unsigned int k = A.pK; k < A.K; k ++)
      Z_r[A.K - k - 1] = reads.uns.v[r].is_compatible_w(negCenters.transcripts[k - A.pK], A);
    uZ.push_back(Z_r);
  }
  for (unsigned int r = 0; r < reads.pos.size; r ++){
    boost::dynamic_bitset<> Z_r(A.K, 0);
    for (unsigned int k = 0; k < A.pK; k ++)
      Z_r[A.K - k - 1] = reads.pos.v[r].is_compatible_w(posCenters.transcripts[k], A);
    pZ.push_back(Z_r);
  }
  for (unsigned int r = 0; r < reads.neg.size; r ++){
    boost::dynamic_bitset<> Z_r(A.K, 0);
    for (unsigned int k = A.pK; k < A.K; k ++)
      Z_r[A.K - k - 1] = reads.neg.v[r].is_compatible_w(negCenters.transcripts[k - A.pK], A);
    nZ.push_back(Z_r);
  }
}

void TranscriptModel::print_Z(vector<boost::dynamic_bitset<> > &aZ){
  for (unsigned int i=0;i<aZ.size();i++)
    cout<<"read "<<i<<": "<< aZ[i]<<"\n";
} 

void TranscriptModel::setX(gsl_matrix* my_X, Reads &reads, vector<boost::dynamic_bitset<> > &sZ, Arguments &A){
  unsigned int z;
  for (unsigned int r = 0; r < reads.size; r ++){
    z = sZ[r].to_ulong();
    for (unsigned int i = 0; i < A.N; i ++)
      gsl_matrix_set(my_X, i, z, gsl_matrix_get(my_X, i, z) + reads.v[r].counts[i]);
  }
}

void TranscriptModel::setX(allReads &reads, Arguments &A){
  setZ(reads,A);
  setX(X, reads.uns, uZ, A);
  setX(X, reads.pos, pZ, A);
  setX(X, reads.neg, nZ, A);
  for (unsigned int i = 0; i < A.N; i ++)
    gsl_matrix_set(X, i, 0, gsl_matrix_get(X, i, 0) + reads.uCounts[i]);
}

void TranscriptModel::setL(Arguments &A){
  for (unsigned int tt = 0; tt < A.pK; tt ++)
    gsl_vector_set(L, tt, posCenters.transcripts[tt].len(A));
  for (unsigned int tt = A.pK; tt < A.K; tt ++)
    gsl_vector_set(L, tt, negCenters.transcripts[tt - A.pK].len(A));
}

//pos: my_AK=A.K
//neg: my_AK=A.K+A.pK
//uns, pos, pZ, pZproposed, 
void TranscriptModel::setZproposed(Reads &uns, Reads &my_pos, vector<boost::dynamic_bitset<> > &my_pZ, vector<boost::dynamic_bitset<> > &my_pZproposed, unsigned int &AK, Arguments &A){
  unsigned int R;
  vector<boost::dynamic_bitset<> > is_modified_u, is_modified_p;
  
  for (unsigned int tt = 0; (unsigned int) tt < proposedTranscripts.size(); tt ++){
    is_modified_u.push_back(boost::dynamic_bitset<> (uns.size, 0));
    is_modified_p.push_back(boost::dynamic_bitset<> (my_pos.size, 0));
  }
  
  uZproposed = uZ;
  my_pZproposed = my_pZ;
  
  for (unsigned int cc = 0; cc < A.lC - 1; cc ++){
    for (unsigned int r = 0; r < uns.InInterval[cc].size(); r ++){
      R = uns.InInterval[cc][r];
      for (unsigned int tt = 0; (unsigned int) tt < proposedTranscripts.size(); tt ++){
	if ( !(is_modified_u[tt][R]) ){//is_modified_u[tt][R]==0 
	  if (proposedTranscripts[tt].modifiedExon[A.lC - cc - 2]){//== 1            
	    is_modified_u.at(tt)[R] = 1;                                                  
	    unsigned int k = AK + proposedTranscripts[tt].t;
	    uZproposed.at(R)[A.K - k - 1] = uns.v[R].is_compatible_w(proposedTranscripts[tt], A);
	  }
	}
      }                                                                     
    }
    //modify Z of positive reads falling in nexon                  
    for (unsigned int r = 0; r < my_pos.InInterval[cc].size(); r ++){
      R = my_pos.InInterval[cc][r];
      for (unsigned int tt = 0; (unsigned int) tt < proposedTranscripts.size(); tt ++){
	  if ( !(is_modified_p[tt][R]) ){//is_modified_s[tt][R]==0   
	    if (proposedTranscripts[tt].modifiedExon[A.lC - cc - 2]){//== 1             
	      is_modified_p.at(tt)[R] = 1;                                       
	      unsigned int k = AK + proposedTranscripts[tt].t;
	      my_pZproposed.at(R)[A.K - k - 1] = my_pos.v[R].is_compatible_w(proposedTranscripts[tt], A);
	    }
	  }                                                           
      }
    }
  }
}

//choose a transcript; choose two centers c and c2; swap indicators of center c and swap indicator of center c2 (for example if center c is included remove it)       
void TranscriptModel::setXproposed(allReads &reads, Reads &my_pos, Reads &my_neg, vector<boost::dynamic_bitset<> > &my_pZ, vector<boost::dynamic_bitset<> > &my_pZproposed, vector<boost::dynamic_bitset<> > &my_nZ, unsigned int AK, Arguments &A){
  setZproposed(reads.uns, my_pos, my_pZ, my_pZproposed, AK, A);
  
  gsl_matrix_set_zero(Xproposed);
  setX(Xproposed, reads.uns, uZproposed, A);
  setX (Xproposed, my_pos, my_pZproposed, A);
  setX (Xproposed, my_neg, my_nZ, A); 
  for (unsigned int i = 0; i < A.N; i ++)
    gsl_matrix_set(Xproposed, i, 0, gsl_matrix_get(Xproposed, i, 0) + reads.uCounts[i]);
}

//AK = 0 or A.pK
void TranscriptModel::setLproposed(unsigned int AK, Arguments &A){//!
  gsl_vector_memcpy(Lproposed, L);
  for (unsigned int p = 0; p < proposedTranscripts.size(); p ++){
    gsl_vector_set(Lproposed, AK + proposedTranscripts[p].t, proposedTranscripts[p].len(A));
  }
}

//set K, posCenters, negCenters, A.pK, A.nK, A.K, A.powK, L, X
TranscriptModel::TranscriptModel(string &genePredFile, allReads &reads, Arguments &A){
  unsigned int a = 0;
  
  if (A.genePrediction==0){
    cerr << "ERROR: option genePrediction = 0 is not implemented"; exit(1);//check this!!!!
  }
  
  while (a == 0){
    a = 1;
    posCenters.init("+", A.pos3, A.pos5, A);        
    negCenters.init("-", A.neg3, A.neg5, A);
    a *= posCenters.check_init(genePredFile, A);
    a *= negCenters.check_init(genePredFile, A);
    if (a == 0) cout<<"Reinitializing Centers\n"; else "Centers initialized\n";
  }
  
  posCenters.init(genePredFile, A);
  negCenters.init(genePredFile, A);
  
  //cout<<"negCenters.K="<<negCenters.K<<"\n";
  //set K, A.K, A.pK, A.nK, A.powK, Z                                            
  A.set(posCenters.K, negCenters.K);
  K = A.K;
  
  //init L and X                       
  L = gsl_vector_alloc(A.K);
  X = gsl_matrix_calloc(A.N, A.powK);
  Lproposed = gsl_vector_alloc(A.K);
  Xproposed = gsl_matrix_alloc(A.N, A.powK);
  
  setL(A);
  setX(reads, A);
  
  posCenters.init_proposal();
  negCenters.init_proposal();
}

TranscriptModel::TranscriptModel(string &genePredFile, Centers &old_posCenters, Centers &old_negCenters, allReads &reads, Arguments &A){
  bool a = 0;
  
  if (A.genePrediction==0){
    cerr << "ERROR: option genePrediction = 0 is not implemented"; exit(1);
  }
  
  posCenters = old_posCenters;
  negCenters = old_negCenters;

  posCenters.clear();
  negCenters.clear();
  
  a *= posCenters.check_init(genePredFile, A);
  a *= negCenters.check_init(genePredFile, A);
  posCenters.init(genePredFile, A);
  negCenters.init(genePredFile, A);
  
  //set K, A.K, A.pK, A.nK, A.powK, Z                         
  A.set2(posCenters.K, negCenters.K);
  K = A.K;
  
  //init L and X                                                                      
  L = gsl_vector_alloc(A.K);
  X = gsl_matrix_calloc(A.N, A.powK);
  //Lproposed = gsl_vector_alloc(A.K);
  //Xproposed = gsl_matrix_alloc(A.N, A.powK);
  
  setL(A);
  setX(reads, A);
}

Move * TranscriptModel::draw_move(ModifiedCenter &modifiedCenter, Arguments &A){//***********************************************
  //draw type of move//center tt    
  //only depends on the center c and on wether or not it is used by transcript t                    
  //if transcript = transcripts[t] and centers[c] = c:    
  //it originally was center.draw_type(transcript.is_using[c], A)                         
  //with probability 10% swap or recombine (if K > 1)
  Center * center = &(Move::centers -> v[modifiedCenter.c]);
  double ProbOrder2t = 0.2, p2t = 1., p_c, u;
  
  if ((Move::centers -> K) > 1)
    p2t = gsl_rng_uniform(A.rnd);
  if (p2t > ProbOrder2t){
    p_c = gsl_ran_beta(A.rnd, 1. + (center -> n_proposed), 1. + (Move::centers -> K) - (center -> n_proposed));
    u = gsl_rng_uniform(A.rnd);
    if (u > p_c){
      if ((center -> is_used_proposed[modifiedCenter.t]) == 0){
	///////cout<<"center -> is_used_proposed[modifiedCenter.t]) == 0\n";
	return new NullMove;
      }else{
	//cout<<"ProposeRemoveCenter\n";
	return new ProposeRemoveCenter;
      }
    }else{
      if ((center -> is_used_proposed[modifiedCenter.t]) == 1){
	if ((center -> ss3Pointers.size) + (center -> ss5Pointers.size) > 2){
	  if (gsl_rng_uniform(A.rnd) < (double) ((center -> ss3Pointers.size) - 1) / (double) ((center -> ss3Pointers.size) + (center -> ss5Pointers.size) - 2)){
	      //cout << "ProposeChangeSS3\n";
	    return new ProposeChangeSS3;
	  }else{
	    //cout << "ProposeChangeSS5\n";
	    return new ProposeChangeSS5;
	  }
	}
      }else{
	//cout<<"ProposeAddCenter\n";
	return new ProposeAddCenter;
      }
      }
    //
  } else if (p2t < 0.05)
	   return new ProposeSwap;
    else
      return new ProposeRecombine;
  return new NullMove;
}

void TranscriptModel::reject(){
  unsigned int c, t;
  
  for (unsigned int i = 0; i < Move::modifiedCenter.size(); i ++){ 
      c = Move::modifiedCenter[i].c;
      t = Move::modifiedCenter[i].t;
      Move::centers -> v.at(c).n_proposed = Move::centers -> v[c].n;
      Move::centers -> v.at(c).v3_proposed.at(t) = Move::centers -> v[c].v3[t];
      Move::centers -> v.at(c).v5_proposed.at(t) = Move::centers -> v[c].v5[t];
      Move::centers -> v.at(c).is_used_proposed = Move::centers -> v[c].is_used;
  }
  for (unsigned int i = 0; i < Move::s.size(); i ++)
    Move::s[i] -> n_proposed = Move::s[i] -> n;
}

void TranscriptModel::accept_reject(allReads &reads, Reads &my_pos, Reads &my_neg, vector<boost::dynamic_bitset<> > &my_pZ, vector<boost::dynamic_bitset<> > &my_pZproposed, vector<boost::dynamic_bitset<> > &my_nZ, unsigned int AK, gsl_vector* lambda, double &epsilon, double &log_epsilon, double &ll, Arguments &A){
  unsigned int c, t;
  double alpha, llproposed;
  
  proposedTranscripts.clear();
  for (unsigned int tt = 0; tt < (Move::t).size(); tt ++){
    t = Move::t[tt];
    ProposedTranscript transcript(t, Move::centers, A);
    transcript.modifiedExon = (transcript.z ^ Move::centers -> transcripts[t].z);
    proposedTranscripts.push_back(transcript);
  }
  //set Lproposed               
  setLproposed(AK, A);
  
  for (unsigned int k = 0; k < A.K; k ++){
    if (gsl_vector_get(Lproposed, k) < 0){
      reject();
      return;
    }
  }
  
  setXproposed(reads, my_pos, my_neg, my_pZ, my_pZproposed, my_nZ, AK, A);    
  llproposed = log_likelihood(lambda, epsilon, log_epsilon, Xproposed, Lproposed, A);
  
  alpha = llproposed;
  alpha -= ll;
  alpha = min(0.,alpha);
  if (gsl_rng_uniform(A.rnd) < exp(alpha)){//accept
    if (Move::small_world == 1)
      Move::ar_counter ++;
    ll = llproposed;
    gsl_vector_swap(L, Lproposed);
    gsl_matrix_swap(X, Xproposed);
    swap(uZ, uZproposed);
    swap(my_pZ, my_pZproposed);
    for (unsigned int i = 0; i < proposedTranscripts.size(); i ++){
      t = proposedTranscripts[i].t;
      swap(Move::centers -> transcripts[t].z, proposedTranscripts[i].z);
    }
    
    for (unsigned int i = 0; i < Move::modifiedCenter.size(); i ++){//copy proposal
      c = Move::modifiedCenter[i].c;
      t = Move::modifiedCenter[i].t;
      Move::centers -> v.at(c).n = Move::centers -> v[c].n_proposed;
      Move::centers -> v.at(c).v3.at(t) = Move::centers -> v[c].v3_proposed[t];
      Move::centers -> v.at(c).v5.at(t) = Move::centers -> v[c].v5_proposed[t];
      Move::centers -> v.at(c).is_used = Move::centers -> v[c].is_used_proposed;
    }
    for (unsigned int i = 0; i < Move::s.size(); i ++)
      Move::s[i] -> n = Move::s[i] -> n_proposed;
    
    //acceptance ratios
    Move::centers -> ar_t ++;
    for (unsigned int i = 0; i < (Move::move_proposed.size()); i ++)
	Move::centers -> ar[Move::move_proposed[i]] += 1;
  }else
    reject();
}

void TranscriptModel::accept_reject_pos(allReads &reads, gsl_vector* lambda, double &epsilon, double &log_epsilon, double &ll, Arguments &A){
  accept_reject(reads, reads.pos, reads.neg, pZ, pZproposed, nZ, 0, lambda, epsilon, log_epsilon, ll, A);
}

void TranscriptModel::accept_reject_neg(allReads &reads, gsl_vector *lambda, double &epsilon, double &log_epsilon, double &ll, Arguments &A){
  accept_reject(reads, reads.neg, reads.pos, nZ, nZproposed, pZ, A.pK, lambda, epsilon, log_epsilon, ll, A);
}

void TranscriptModel::propose_update(unsigned int &t1, allReads &reads, gsl_vector *lambda, double &epsilon, double &log_epsilon, double &ll, unsigned int &K, Arguments &A){
  Move::clearMove();
  Move::small_world = 0;
  //////////////////// draw t
  unsigned int c1;
  double u = gsl_rng_uniform(A.rnd);
  (Move::centers -> ar_den_t) += 1;
  ////////////////// draw c
  c1 = gsl_rng_uniform_int(A.rnd, Move::centers -> size);
  ////////////////// propose move
  ModifiedCenter modifiedCenter(c1, t1);
  Move *move = draw_move(modifiedCenter, A);
  move -> propose(modifiedCenter, A);
  delete move;
  /////////////////////////////////////////////////////////////////      
  if (K > 1){
    if (u > 0.7){
      unsigned int t2 = gsl_rng_uniform_int(A.rnd, K - 1);
      if (t2 == t1) 
	t2 = K - 1;
      ModifiedCenter modifiedCenter2(c1, t2);
      Move *move2 = draw_move(modifiedCenter2, A);
      move2 -> propose(modifiedCenter2, A);
      delete move2;
    }
  }
  
  double u_sw = gsl_rng_uniform(A.rnd);
  if (u_sw > 0.7){////////////////0.9
    Move::small_world = 1;
    Move::counter ++;
    unsigned int c2;
    for (unsigned int i = 0; i < 5; i ++){
      c2 = gsl_rng_uniform_int(A.rnd, Move::centers -> size);
      ModifiedCenter modifiedCenter3(c2, t1);
      Move *move3 = draw_move(modifiedCenter3, A);
      move3 -> propose(modifiedCenter3, A);
      delete move3;
    }
  }
  Move::sortMove();
}

void TranscriptModel::update(allReads &reads, gsl_vector *lambda, double &epsilon, double &log_epsilon, double &ll, Arguments &A){
  unsigned int t1 = gsl_rng_uniform_int(A.rnd, A.K);
  
  if (t1 < A.pK ){
    Move::centers = &posCenters;
    propose_update(t1, reads, lambda, epsilon, log_epsilon, ll, A.pK, A);
    if ((Move::modifiedCenter.size()) != 0)
      accept_reject_pos(reads, lambda, epsilon, log_epsilon, ll, A);
  }else{
    t1 -= A.pK;
    Move::centers = &negCenters;
    propose_update(t1, reads, lambda, epsilon, log_epsilon, ll, A.nK, A);
    if ((Move::modifiedCenter.size()) != 0)
      accept_reject_neg(reads, lambda, epsilon, log_epsilon, ll, A);
  }
}

TranscriptModel::~TranscriptModel(){
  gsl_matrix_free(Xproposed);
  gsl_vector_free(Lproposed);
  gsl_matrix_free(X);
  gsl_vector_free(L);
}

void TranscriptModel::print_X(){
  print_my_gsl_matrix(X); 
}

void TranscriptModel::print_L(){
  print_my_gsl_vector(L);
}

void TranscriptModel::print(string strand, Centers &centers, Arguments &A){
  if (centers.K > 0){
    //I need to make a copy of the transcripts because I'm changing their order
    ////////////////////////////////////////
    vector<Transcript> my_transcripts = centers.transcripts;
    vector<Transcript>::iterator it;
    A.buffer << strand << ",";
    sort(my_transcripts.begin(), my_transcripts.end());
    it = unique(my_transcripts.begin(), my_transcripts.end());
    my_transcripts.resize( it - my_transcripts.begin() ); 
    
    for (it = my_transcripts.begin(); it != my_transcripts.end(); ++ it)
      A.buffer << it -> z << ",";
      A.buffer << " ";
  }
}

void TranscriptModel::print_to_Abuffer(const double &ll, Arguments &A){
  print("+", posCenters, A);
  print("-", negCenters, A);
  A.buffer << fixed << ll << "\n";
}

