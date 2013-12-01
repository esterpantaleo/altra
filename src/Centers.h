/** \file Centers.h
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

#ifndef CENTERS_H
#define CENTERS_H

#include "Transcript.h"

class Centers{
public:
  unsigned int size, K;
  string strand;
  vector<SS> three_primes, five_primes;
  vector<Center> v;// contains the set of centers, of size "size"          
  vector<Transcript> transcripts;
  double ar_t, ar_den_t;
  vector<double> ar, ar_den; 

  void do_nothing(){}
  void clear();
  void print_ar(string strand, Arguments &A);
  void clear_ar(Arguments &A);
  void init_proposal();
  void init(string my_strand, vector<unsigned int> &pos3, vector<unsigned int> &pos5, Arguments &A);
  bool get_c3_c5(unsigned int &genePred5SS, unsigned int &genePred3SS, unsigned int &c3, unsigned int &c5, Arguments &A);
  void init2(vector<unsigned int> &genePred3SS, vector<unsigned int> &genePred5SS, Arguments &A);
  bool check_init(string &genePredFile, Arguments &A);
  void init(string &genePredFile, Arguments &A);
};

#endif
