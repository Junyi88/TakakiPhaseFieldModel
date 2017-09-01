#ifndef TAKPHASE_H
#define TAKPHASE_H

#include "JFDMpi.h"

template <class FDClass>
class ACBulk {
public:
  ACBulk(JMpi inJMpi, const JMat &in_F);
  ACBulk(JMpi inJMpi);
	ACBulk & operator= (const ACBulk &in1); //Write to operator

  //=========================
  // Useful Functions
  // void SetF(const JMat &in_F);
  // void Calc_All();
  // void Calc_FD();
  // void Calc_Powers();

protected:
  JMpi _MpiObj;
	int _NY, _NX, _Ny;

  JMat _F;
  JMat _F2, _F3, _F4;
  FDClass _D;

  MPI_Status _status;
};


#endif
