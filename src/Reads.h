/** \file Reads.h
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

#ifndef READS_H
#define READS_H

#include "Transcript.h"

class Read{
public:
  unsigned int size;
  vector<unsigned int> c5, c3; //5' ends of read, 3' ends of read
  vector <double> counts;

  Read(Arguments &A);
  Read(string &my_string, Arguments &A);
  bool is_compatible_w(Transcript &transcript, Arguments &A);
};

class Reads{
public:
  int size;
  vector<vector<int> > InInterval;
  vector<Read> v;
  
  void set_zero(Arguments &A);
  void get(string &my_string, int &individual, Arguments &A);
  void print();
  void print_InInterval();
};


//get the reads from the shell: pReads=positive reads (spliced reads containing a junction with a positive sense splice signal); 
//nReads=negative reads; uRead=unsigned reads (unspliced reads, that map to the reference sequence)
class allReads{
public:
  Reads pos;
  Reads neg;
  Reads uns;
  vector<int> uCounts; // number of reads that cannot map to any transcript  

  allReads (Arguments &A);
};

#endif
