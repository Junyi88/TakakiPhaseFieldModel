#ifndef GClass_H
#define GClass_H

#include "JMat.h"

template <class T>
class GClass{
public:
  GClass(const T &in_F) : F(in_F){};
  T _F;
};

#endif
