#ifndef TAKACWALLENERGY_H
#define TAKACWALLENERGY_H

#include "TakPhase.h"

//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
class TakACWallEnergy {
public:
  TakACWallEnergy(TakPhase<FDClass> * inPhi, const double &inWa, JMpi inJMpi);
	TakACWallEnergy<FDClass> & operator= (const TakACWallEnergy<FDClass> &in1); //Write to operator

  void Calc_All();

  // Getter Functions
  TakPhase PhiP(){return _Phi;};
  double dFdPhase(const int &y, const int &x) {return _dFdphase(y,x);};
  JMat * dFdPhasePointer() {return &(_dFdPhase);};

protected:
  TakPhase<FDClass> * _Phi;
  JMpi _MpiObj;
	int _NY, _NX, _Ny;

  double _Wa;
  JMat _dFdPhase;
};

#endif
