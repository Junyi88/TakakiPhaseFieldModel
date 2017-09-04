#ifndef GClass_H
#define GClass_H

#include "JMat.h"

template <class T>
class GClass{
public:
  // GClass(const T &in_F) : _F(in_F){};
  GClass(const T &in_F){};
  T _F;
};


GClass<class T>::GClass(const T &in_F) : _F(in_F){};
#endif
