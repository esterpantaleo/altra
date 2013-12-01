/** \file Move.h
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

#ifndef MOVE_H
#define MOVE_H

#include "Centers.h"

extern enum TypeOfMove{
  PChangeSS5,
  PChangeSS3,
  PAddCenter,
  PRemoveCenter,
  PSwap,
  PRecombine,
  PNull
} Types;

//if element a and b are the same return 1                                                   
template<class T> bool swap(vector<T> &v, const unsigned int &a, const unsigned int &b);

class Move{
public:
  static vector<TypeOfMove> move_proposed;
  static vector<ModifiedCenter> modifiedCenter;
  static vector<unsigned int> t;
  static vector<SS*> s;
  static Centers* centers;
  static int counter;
  static int ar_counter; //small_world acceptance ratio
  static int small_world;

  virtual void propose(ModifiedCenter &modCenter, Arguments &A) = 0;
  static void clearMove();
  static void sortMove();
  void proposeChangeSS(ModifiedCenter &modCenter, SS* &my_ss, SSpointers & ss, Arguments &A);
  void proposeSR(unsigned int &t1, unsigned int &c1, unsigned int c2, Arguments &A);
};

class NullMove: public Move{
public:
  void propose(ModifiedCenter &modCenter, Arguments &A);
};

class ProposeChangeSS5: public Move{
public:
  void propose(ModifiedCenter &modCenter, Arguments &A);
};

class ProposeChangeSS3: public Move{
public:
  void propose(ModifiedCenter &modCenter, Arguments &A);
};

class ProposeAddCenter: public Move{
public:
  void propose(ModifiedCenter &modCenter, Arguments &A);
};

class ProposeRemoveCenter: public Move{
public:
  void propose(ModifiedCenter &modCenter, Arguments &A);
};

class ProposeRecombine: public Move{
public:
  void propose(ModifiedCenter &modCenter, Arguments &A);
};

class ProposeSwap: public Move{
public:
  void propose(ModifiedCenter &modCenter, Arguments &A);
};
 
#endif
