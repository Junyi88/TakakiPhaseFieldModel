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
	int _NY, _NX, _Ny;

  JMat _F;
  JMat _F2, _F3, _F4, _F5;
  FDClass _D;
  JMat _P;
  JMat _dP;

  MPI_Status _status;
};


#endif
