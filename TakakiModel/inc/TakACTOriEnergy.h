#ifndef TAKACTORIENERGY_H
#define TAKACTORIENERGY_H

#include "TakPhase.h"
#include "TakAngle.h"

//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass, class FDAngleClass>
class TakACTOriEnergy {
public:
  TakACTOriEnergy(TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta,
    const double &insConst, const double &inMTheta0,  const double &inInvPhiMin,
    JMpi inJMpi);
	TakACTOriEnergy<FDClass, FDAngleClass> & operator= (const TakACTOriEnergy<FDClass, FDAngleClass> &in1); //Write to operator

  void Calc_dFdPhase();

  void Calc_MTheta();
  void Calc_dAngleFront();
  void Calc_dAngleRear();
  void Calc_dThetadt();

  void Calc_All();

  // Getter Functions

  double dFdPhase(const int &y, const int &x) {return _dFdPhase(y,x);};
  double dThetadt(const int &y, const int &x) {return _dThetadt(y,x);};

protected:
  TakPhase<FDClass> * _Phi;
  TakAngle<FDAngleClass> * _Theta;
  JMpi _MpiObj;
	int _NY, _NX, _Ny;

  double _sConst, _MTheta0, _Ms, _s2, _InvPhiMin;
  JMat _dFdPhase;

  JMat _MTheta;
  JMat _dAngleFront;
  JMat _dAngleRear;
  JMat _dThetadt;
};

#endif
