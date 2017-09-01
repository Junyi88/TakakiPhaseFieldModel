#ifndef TAKPHASE_H
#define TAKPHASE_H

#include "JFDMpi.h"

template <class FDClass>
class TakPhase {
public:
  TakPhase(JMpi inJMpi, const JMat &in_F);
  TakPhase(JMpi inJMpi);
	TakPhase<FDClass> & operator= (const TakPhase<FDClass> &in1); //Write to operator

  //=========================
  // Useful Functions
  // void SetF(const JMat &in_F); // DANGER DON'T USE

  void Calc_FD();
  void Calc_Powers();
  void Calc_P();
  void Calc_dP();
  void Calc_All();

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
