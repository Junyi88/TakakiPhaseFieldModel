#ifndef TAKPHASE_H
#define TAKPHASE_H

#include "JFDMpi.h"

//## ===================================
template <class FDClass>
class TakPhase {
public:
  TakPhase(JMpi inJMpi, const JMat &in_F);

protected:
  JMpi _MpiObj;


  FDClass _F;
  
};


#endif
