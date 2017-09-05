#include "TakACWallEnergy.h"

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
