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

//***********************************************************************
//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
TakACGradEnergy<FDClass>::TakACGradEnergy(TakPhase<FDClass> * inPhi, const double &inalpha, JMpi inJMpi) :
 _Phi(inPhi), _MpiObj(inJMpi), _NY(_MpiObj.NYGl()),  _NX(_MpiObj.NX()),
 _Ny(_MpiObj.NYLo()), _alpha(inalpha), _alpha2(2.0*_alpha), _dFdPhase(_NY,_NX)  {}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
TakACGradEnergy<FDClass> & TakACGradEnergy<FDClass>::operator= (const TakACGradEnergy<FDClass> &in1) {
   _Phi=in1.PhiP();
   _MpiObj=in1._MpiObj;
   _NY=in1._NY;
   _NX=in1._NX;
   _Ny=in1._Ny;
   _alpha=in1._alpha;
   _alpha2=in1._alpha2;
   _dFdPhase=in1._dFdPhase;
 }

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
void TakACGradEnergy<FDClass>::Calc_All(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _dFdPhase(j,i)=_Phi->D2(j,i);
      _dFdPhase(j,i)*=-_alpha2;
    }
}

#endif
