#ifndef TAKACGRADENERGY_H
#define TAKACGRADENERGY_H

#include "TakPhase.h"

//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
class TakACGradEnergy {
public:
  TakACGradEnergy(TakPhase<FDClass> * inPhi, const double &inalpha, JMpi inJMpi);
	TakACGradEnergy<FDClass> & operator= (const TakACGradEnergy<FDClass> &in1); //Write to operator

  void Calc_All();

  // Getter Functions
  TakPhase<FDClass> * PhiP(){return &(_Phi);};
  double dFdPhase(const int &y, const int &x) {return _dFdPhase(y,x);};
  JMat * dFdPhasePointer() {return &(_dFdPhase);};

protected:
  TakPhase<FDClass> * _Phi;
  JMpi _MpiObj;
	int _NY, _NX, _Ny;

  double _alpha;
  double _alpha2;
  JMat _dFdPhase;
};

#endif
