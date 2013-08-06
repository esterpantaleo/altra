/** \file Center.cpp
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

#include "Center.h"

SS::SS(unsigned int &cc){
  c = cc;
  n = 0;
  n_proposed = 0;                                                         
}

bool SS_less(const SS* lhs, const SS* rhs){
  return(lhs -> c < rhs -> c);
}

bool SS_equal(const SS* lhs, const SS* rhs){
  return(lhs -> c  == rhs -> c);
}

void Center::clear(){
  for (unsigned int i = 0; i < ss3Pointers.size; i ++)
    ss3Pointers.v[i] -> n = 0;
  for (unsigned int i = 0; i < ss5Pointers.size; i ++)
    ss5Pointers.v[i] -> n = 0;
  v5.clear();
  v3.clear();
  n = 0;
  is_used.clear();
}

//define ss3Pointers.v and ss3Pointers.size 
Center::Center(unsigned int &i3, unsigned int &j5, vector<SS> &ss3, vector<SS> &ss5, Arguments &A){
  //define a center between A.listCoordinate[ss3[i3].c] and A.listCoordinate[ss5[j5].c]
  //define center's possible 3' and 5' ss (in the vector v), and their sizes
  vector<SS*>::iterator it;
  
  for (unsigned int i = 0; i <= i3; i ++){
    for (unsigned int j = j5; j < ss5.size(); j ++){
      if (A.listCoordinate[ss5[j].c] - A.listCoordinate[ss3[i3 - i].c] <= 2 * A.MAX_EX_LEN_DIV_2){
	if (A.listCoordinate[ss5[j].c] - A.listCoordinate[ss3[i3 - i].c] > A.MIN_EX_LEN){
	  ss3Pointers.v.push_back(&(ss3[i3 - i]));
	  ss5Pointers.v.push_back(&(ss5[j]));
	}
      }else
	break;
    }
  }
  
  sort(ss3Pointers.v.begin(), ss3Pointers.v.end(), SS_less);
  it = unique(ss3Pointers.v.begin(), ss3Pointers.v.end(), SS_equal);
  ss3Pointers.v.resize( it - ss3Pointers.v.begin() );
  sort(ss5Pointers.v.begin(), ss5Pointers.v.end(), SS_less);
  it = unique(ss5Pointers.v.begin(), ss5Pointers.v.end(), SS_equal);
  ss5Pointers.v.resize(it - ss5Pointers.v.begin());
  ss3Pointers.size = (unsigned int) ss3Pointers.v.size();
  ss5Pointers.size = (unsigned int) ss5Pointers.v.size();
  n = 0;
}

//get the indices i3 (and i5) of the  elements in three_primes (five_primes) that have value the_3 (the_5) 
//return 1(0) if exon is (not) compatible with center
bool Center::get(unsigned int &c3, unsigned int &c5, SS*& s3, SS*& s5){
  unsigned int counter = 0, p3, p5;
  for (p3 = 0; p3 < ss3Pointers.size; p3 ++){
    s3 = ss3Pointers.v[p3];
    if ((s3 -> c) == c3){
      counter = 1;
      break;
    }
  }
  for (p5 = 0; p5 < ss5Pointers.size; p5 ++){
    s5 = ss5Pointers.v[p5];
    if ((s5 -> c) == c5){
      counter += 1;
      break;
    }
  }
  if (counter == 2)
    return 1;
  return 0;
}

SS* Center::draw_SS_from_posterior(SSpointers &ssPointers, Arguments &A){
  //possible_3=3,4,7,8                                                      
  //unsigned int my_proposed_ss; //chosen  SS                                       
  double q[ssPointers.size], alpha[ssPointers.size], q_cumul, u = gsl_rng_uniform(A.rnd);
  //sample q form the posterior (dirichlet)                                 
  for (unsigned int s = 0; s < ssPointers.size; s ++)
    alpha[s] = A.alphaSS + (ssPointers.v[s] -> n_proposed);
  gsl_ran_dirichlet(A.rnd, ssPointers.size, alpha, q);
  //propose SS from the posterior                                 
  q_cumul = q[0];
  for (unsigned int s = 0; s < ssPointers.size; s ++){
    if (u < q_cumul)
      return ssPointers.v[s];
    q_cumul += q[s + 1];
  }
}

void Center::init_not_used(Arguments &A){
  int u = gsl_rng_uniform_int(A.rnd, ss5Pointers.size), uu = gsl_rng_uniform_int(A.rnd, ss3Pointers.size);
  is_used.push_back(0);
  v5.push_back(ss5Pointers.v[u]);
  v3.push_back(ss3Pointers.v[uu]);
}

ModifiedCenter::ModifiedCenter(unsigned int &my_c, unsigned int &my_t){
  t = my_t;
  c = my_c;
}

bool is_less_than(ModifiedCenter x, ModifiedCenter y){
  return ((x.t == y.t && x.c < y.c) || x.t < y.t);
}

bool is_equal_to(ModifiedCenter x, ModifiedCenter y){
  return ((x.c == y.c) && (x.t == y.t));
}

