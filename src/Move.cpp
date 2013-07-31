#include "Move.h"

template<class T> bool swap(vector<T> &v, const unsigned int &a, const unsigned int &b){
  T tmp = v[a];

  if (tmp == v[b])
    return 1;
  v.at(a) = v[b];
  v.at(b) = tmp;

  return 0;
}

void Move::clearMove(){
  move_proposed.clear();
  modifiedCenter.clear();
  t.clear();
  s.clear();
}

void Move::sortMove(){
  vector<SS*>::iterator it;
  vector<unsigned int>::iterator It;
  vector<ModifiedCenter>::iterator iT;
  
  sort(s.begin(), s.end(), SS_less);
  it = unique(s.begin(), s.end(), SS_equal);
  s.resize(it - s.begin());  
  sort(t.begin(), t.end());
  It = unique(t.begin(), t.end());
  t.resize(It - t.begin());
  sort(modifiedCenter.begin(), modifiedCenter.end(), is_less_than);
  iT = unique(modifiedCenter.begin(), modifiedCenter.end(), is_equal_to);
  modifiedCenter.resize(iT - modifiedCenter.begin());
}

void Move::proposeChangeSS(ModifiedCenter &modCenter, SS* &my_ss, SSpointers & ss, Arguments &A){
  SS* s_current = my_ss;
  
  my_ss = centers -> v.at(modCenter.c).draw_SS_from_posterior(ss, A);
  if (s_current != my_ss){//not accepted yet 
    (my_ss -> n_proposed) += 1;
    (s_current -> n_proposed) -= 1;
    
    s.push_back(my_ss);
    s.push_back(s_current);
    t.push_back(modCenter.t);
    modifiedCenter.push_back(modCenter);
  }
}

void Move::proposeSR(unsigned int &t1, unsigned int &c1, unsigned int c2, Arguments &A){
  //if proposing a recombination of all centers you are proposing to just swap the two transcripts 
  //which is essentially a null move
  if (c1 - c2 == (Move::centers -> size))
    return;
  
  bool are_same, is_same;
  //sample a transcript t2 uniformly in [0,K-1] where K is the number of transcripts
  //if t2==t1 then take the K-th transcript
  unsigned int one_different = 0, t2 = gsl_rng_uniform_int(A.rnd, (centers -> K) - 1);
  
  if (t2 == t1)
    t2 = (centers -> K) - 1;
  for (unsigned int cc = c2; cc <= c1; cc ++){//if one of them is different (is_same==0) then they are different (are_same=0)
    are_same = 1;
    is_same = swap(centers -> v.at(cc).v3_proposed, t1, t2);
    are_same *= is_same;
    is_same = swap(centers -> v.at(cc).v5_proposed, t1, t2);
    are_same *= is_same; 
    is_same = swap(centers -> v.at(cc).is_used_proposed, t1, t2);
    are_same *= is_same;
    if (are_same == 0){
      one_different ++;
      modifiedCenter.push_back(ModifiedCenter(cc, t1));
      modifiedCenter.push_back(ModifiedCenter(cc, t2));
    }
  }
  if (one_different != 0){
    t.push_back(t1);
    t.push_back(t2);
  }
}

void NullMove::propose(ModifiedCenter &modCenter, Arguments &A){
  move_proposed.push_back(PNull);
}

void ProposeChangeSS5::propose(ModifiedCenter &modCenter, Arguments &A){
  proposeChangeSS(modCenter, centers -> v.at(modCenter.c).v5_proposed.at(modCenter.t), centers -> v.at(modCenter.c).ss5Pointers, A);
  centers -> ar_den[PChangeSS5] += 1;
  move_proposed.push_back(PChangeSS5);                    
}

void ProposeChangeSS3::propose(ModifiedCenter &modCenter, Arguments &A){
    proposeChangeSS(modCenter, centers -> v.at(modCenter.c).v3_proposed.at(modCenter.t), centers -> v.at(modCenter.c).ss3Pointers, A);
    centers -> ar_den[PChangeSS3] += 1;
    move_proposed.push_back(PChangeSS3); 
}

void ProposeAddCenter::propose(ModifiedCenter &modCenter, Arguments &A){
  SS *s5 = centers -> v.at(modCenter.c).draw_SS_from_posterior(centers -> v[modCenter.c].ss5Pointers, A);
  SS *s3 = centers -> v.at(modCenter.c).draw_SS_from_posterior(centers -> v[modCenter.c].ss3Pointers, A);
  s3 -> n_proposed += 1;
  s5 -> n_proposed += 1;
  
  centers -> v[modCenter.c].n_proposed += 1;
  centers -> v[modCenter.c].is_used_proposed.at(modCenter.t) = 1;
  centers -> v.at(modCenter.c).v3_proposed.at(modCenter.t) = s3;
  centers -> v.at(modCenter.c).v5_proposed.at(modCenter.t) = s5;
  s.push_back(s3);
  s.push_back(s5);
  t.push_back(modCenter.t);
  modifiedCenter.push_back(modCenter);
  centers -> ar_den[PAddCenter] += 1;
  move_proposed.push_back(PAddCenter);
}

void ProposeRemoveCenter::propose(ModifiedCenter &modCenter, Arguments &A){
  centers -> v[modCenter.c].n_proposed -= 1;
  centers -> v[modCenter.c].is_used_proposed.at(modCenter.t) = 0;
  centers -> v[modCenter.c].v3_proposed[modCenter.t] -> n_proposed -= 1;
  centers -> v[modCenter.c].v5_proposed[modCenter.t] -> n_proposed -= 1;
  s.push_back(centers -> v[modCenter.c].v3_proposed[modCenter.t]);
  s.push_back(centers -> v[modCenter.c].v5_proposed[modCenter.t]);
  t.push_back(modCenter.t);
  modifiedCenter.push_back(modCenter);
    centers -> ar_den[PRemoveCenter] += 1;
    move_proposed.push_back(PRemoveCenter);
}

void ProposeRecombine::propose(ModifiedCenter &modCenter, Arguments &A){
  proposeSR(modCenter.t, modCenter.c, 0, A);
  centers -> ar_den[PRecombine] += 1;
  move_proposed.push_back(PRecombine);
}

void ProposeSwap::propose(ModifiedCenter &modCenter, Arguments &A){
  proposeSR(modCenter.t, modCenter.c, modCenter.c, A);
  centers -> ar_den[PSwap] += 1;
  move_proposed.push_back(PSwap);                                                  
} 
 
