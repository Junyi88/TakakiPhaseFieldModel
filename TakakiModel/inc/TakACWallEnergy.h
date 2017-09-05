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
  TakPhase<FDClass> * PhiP(){return &(_Phi);};
  double dFdPhase(const int &y, const int &x) {return _dFdPhase(y,x);};
  JMat * dFdPhasePointer() {return &(_dFdPhase);};

protected:
  TakPhase<FDClass> * _Phi;
  JMpi _MpiObj;
	int _NY, _NX, _Ny;

  double _Wa;
  JMat _dFdPhase;
};

//******************************************************************
//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
TakACWallEnergy<FDClass>::TakACWallEnergy(TakPhase<FDClass> * inPhi, const double &inWa, JMpi inJMpi) :
 _Phi(inPhi), _MpiObj(inJMpi), _NY(_MpiObj.NYGl()),  _NX(_MpiObj.NX()),
 _Ny(_MpiObj.NYLo()), _Wa(inWa), _dFdPhase(_NY,_NX)  {}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
TakACWallEnergy<FDClass> & TakACWallEnergy<FDClass>::operator= (const TakACWallEnergy<FDClass> &in1) {
   _Phi=in1.PhiP();
   _MpiObj=in1._MpiObj;
   _NY=in1._NY;
   _NX=in1._NX;
   _Ny=in1._Ny;
   _Wa=in1._Wa;
   _dFdPhase=in1._dFdPhase;
 }

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
void TakACWallEnergy<FDClass>::Calc_All(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _dFdPhase(j,i)=(4.0*(_Phi->F3(j,i))) - (6.0*(_Phi->F2(j,i))) +(2.0*(_Phi->F(j,i)));
      _dFdPhase(j,i)*=_Wa;
    }
}

#endif
