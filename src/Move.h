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
