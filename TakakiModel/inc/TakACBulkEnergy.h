#ifndef TAKACBULKENERGY_H
#define TAKACBULKENERGY_H

#include "TakPhase.h"

//##========================================================================
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

//****************************************************************************
#include "TakACBulkEnergy.h"

//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
TakACBulkEnergy<FDClass>::TakACBulkEnergy(TakPhase<FDClass> * inPhi, const double &BVec,
  const double &mu, const JMat &RhoIn, JMpi inJMpi) :
 _Phi(inPhi), _MpiObj(inJMpi), _NY(_MpiObj.NYGl()),  _NX(_MpiObj.NX()),
 _Ny(_MpiObj.NYLo()), _BVector(BVec), _mu(mu), _Coeff(0.5*_BVector*_BVector*_mu),
 _rho(RhoIn), _EStored(_NY,_NX), _dFdPhase(_NY,_NX)  {
   Calc_EStored();
 }

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
TakACBulkEnergy<FDClass> & TakACBulkEnergy<FDClass>::operator= (const TakACBulkEnergy<FDClass> &in1) {
   _Phi=in1.PhiP();
   _MpiObj=in1._MpiObj;
   _NY=in1._NY;
   _NX=in1._NX;
   _Ny=in1._Ny;

  _BVector=in1._BVector;
  _mu=in1._mu;
  _Coeff=in1._Coeff;

   _rho=in1._rho;
   _EStored=in1._EStored;
   _dFdPhase=in1._dFdPhase;
 }

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
void TakACBulkEnergy<FDClass>::Calc_All(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _dFdPhase(j,i)=-_EStored(j,i)*(_Phi->dP(j,i));
    }
}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
void TakACBulkEnergy<FDClass>::Calc_EStored(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _EStored(j,i)=_Coeff*_rho(j,i);
    }
}


#endif
