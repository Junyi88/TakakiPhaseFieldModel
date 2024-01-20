#ifndef CHENYUNACBULKENERGY_H
#define CHENYUNACBULKENERGY_H

#include "TakPhase.h"

//##========================================================================
template <class FDClass>
class ChenYunACBulkEnergy {
public:
  ChenYunACBulkEnergy(TakPhase<FDClass> * inPhi, const double &epsilon, const double &a, JMpi inJMpi);
	ChenYunACBulkEnergy<FDClass> & operator= (const ChenYunACBulkEnergy<FDClass> &in1); //Write to operator

  void Calc_All();

  // Getter Functions
  TakPhase<FDClass> * PhiP(){return &(_Phi);};
  double dFdPhase(const int &y, const int &x) {return _dFdPhase(y,x);};
  JMat * dFdPhasePointer() {return &(_dFdPhase);};

protected:
  TakPhase<FDClass> * _Phi;
  JMpi _MpiObj;
	int _NY, _NX, _Ny;

  double _epsilon2;
  double _a2;
  JMat _dFdPhase;
};

//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
ChenYunACBulkEnergy<FDClass>::ChenYunACBulkEnergy(TakPhase<FDClass> * inPhi,
  const double& epsilon, const double &a, JMpi inJMpi) :
 _Phi(inPhi), _MpiObj(inJMpi), _NY(_MpiObj.NYGl()),  _NX(_MpiObj.NX()),
 _Ny(_MpiObj.NYLo()), _epsilon2(epsilon), _a2(a * a), _dFdPhase(_NY,_NX)  {
 }

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
ChenYunACBulkEnergy<FDClass> & ChenYunACBulkEnergy<FDClass>::operator= (const ChenYunACBulkEnergy<FDClass> &in1) {
   _Phi=in1.PhiP();
   _MpiObj=in1._MpiObj;
   _NY=in1._NY;
   _NX=in1._NX;
   _Ny=in1._Ny;

   _epsilon2=in1._epsilon2;

   _dFdPhase=in1._dFdPhase;

   return *this;
 }

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
void ChenYunACBulkEnergy<FDClass>::Calc_All(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _dFdPhase(j,i) = _epsilon2 * _Phi->D2(j,i) + (1.0 - _Phi->F(j,i)) * _a2;
    }
}


#endif
