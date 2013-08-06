/** \file Transcript.h
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

#ifndef TRANSCRIPT_H
#define TRANSCRIPT_H

#include "Center.h"

class Transcript{
public:
  unsigned int t;
  boost::dynamic_bitset<> z;
  boost::dynamic_bitset<> modifiedExon;

  bool operator< (const Transcript &rhs) const;
  bool operator== (const Transcript &rhs) const;
  virtual void to_z(vector<Center> &v, unsigned int &size, Arguments &A);

  Transcript(){}  
  Transcript(unsigned int &my_t, Arguments &A);
  Transcript(Arguments &A);
  Transcript(unsigned int &my_t, vector<Center> &v, unsigned int &size, Arguments &A);
  int len(Arguments &A);
};

#endif
