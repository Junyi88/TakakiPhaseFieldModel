#ifndef TAKACWALLENERGY_H
#define TAKACWALLENERGY_H

#include "JFPhase.h"

//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
class TakACWallEnergy {
public:
  TakACWallEnergy(TakPhase<FDClass> * inPhi, const double &inWa, JMpi inJMpi);
	TakACWallEnergy<FDClass> & operator= (const TakACWallEnergy<FDClass> &in1); //Write to operator

  void Calc_All();

  TakPhase PhiP(){return _Phi;};

protected:
  TakPhase<FDClass> * _Phi;
  JMpi _MpiObj;
	int _NY, _NX, _Ny;

  double _Wa;
  JMat _dFdPhase;
};

#endif
