#ifndef TAKACBULKENERGY_H
#define TAKACBULKENERGY_H

#include "TakPhase.h"

//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
class TakACBulkEnergy {
public:
  TakACBulkEnergy(TakPhase<FDClass> * inPhi, const double &BVec, const double &mu,
    const JMat &RhoIn, JMpi inJMpi);
	TakACBulkEnergy<FDClass> & operator= (const TakACBulkEnergy<FDClass> &in1); //Write to operator

  void Calc_EStored();
  void Calc_All();

  // Getter Functions
  TakPhase<FDClass> * PhiP(){return &(_Phi);};
  double dFdPhase(const int &y, const int &x) {return _dFdPhase(y,x);};
  JMat * dFdPhasePointer() {return &(_dFdPhase);};

protected:
  TakPhase<FDClass> * _Phi;
  JMpi _MpiObj;
	int _NY, _NX, _Ny;

  double _BVector, _mu, _Coeff;
  JMat _rho;
  JMat _EStored;
  JMat _dFdPhase;
};

#endif
