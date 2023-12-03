#ifndef ChenYunACTOriEnergy_H
#define ChenYunACTOriEnergy_H

#include "TakPhase.h"
#include "TakAngle.h"

//##========================================================================
template <class FDClass, class FDAngleClass>
class ChenYunACTOriEnergy {
public:
  ChenYunACTOriEnergy(TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta,
    const double& alpha, const double& beta, const double& omega, 
    const double& mu,
    JMpi inJMpi);
	ChenYunACTOriEnergy<FDClass, FDAngleClass> & operator= (const ChenYunACTOriEnergy<FDClass, FDAngleClass> &in1); //Write to operator


  void Calc_Simples();
  void Calc_dThetadt_T1();
  void Calc_dThetadt_T2();
  void Calc_dThetadt();

  void Calc_dFdPhase();
  void Calc_All();

  // Getter Functions

  double dFdPhase(const int &y, const int &x) {return _dFdPhase(y,x);};
  double dThetadt(const int &y, const int &x) {return _dThetadt(y,x);};

  JMat * dFdPhasePointer() {return &(_dFdPhase);};
  JMat * dThetadtPointer() {return &(_dThetadt);};

protected:
  TakPhase<FDClass> * _Phi;
  TakAngle<FDAngleClass> * _Theta;
  JMpi _MpiObj;
	int _NY, _NX, _Ny;

  double _alpha, _beta, _omega, _omega2, _mu;

  JMat _Q; 
  
  JMat _dFdPhase_T1;
  JMat _dFdPhase_T2;
  JMat _dFdPhase;

  JMat _dThetadt_T1;
  JMat _dThetadt_T2;
  JMat _dThetadt;
};

//**********************************************************************
//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass, class FDAngleClass>
ChenYunACTOriEnergy<FDClass, FDAngleClass>::ChenYunACTOriEnergy(
  TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta,
  const double& alpha, const double& beta, const double& omega, 
  const double& mu,
  JMpi inJMpi) : _Phi(inPhi), _Theta(inTheta), _MpiObj(inJMpi),
  _NY(_MpiObj.NYGl()),  _NX(_MpiObj.NX()),  _Ny(_MpiObj.NYLo()),
  _alpha(alpha), _beta(beta), 
  _omega(omega), _omega2(omega*omega), _mu(mu),
  _Q(_NY,_NX), 
  _dFdPhase_T1(_NY,_NX), _dFdPhase_T2(_NY,_NX), _dFdPhase(_NY,_NX),
  _dThetadt_T1(_NY,_NX), _dThetadt_T2(_NY,_NX), _dThetadt(_NY,_NX)
  {}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass, class FDAngleClass>
ChenYunACTOriEnergy<FDClass, FDAngleClass> & ChenYunACTOriEnergy<FDClass, FDAngleClass>::operator=
  (const ChenYunACTOriEnergy<FDClass, FDAngleClass> &in1) {

    _Phi=in1._Phi;
    _Theta=in1._Theta;
    _MpiObj=in1._MpiObj;
  	_NY=in1._NY;
    _NX=in1._NX;
    _Ny=in1._Ny;

    _alpha=in1._alpha;
    _beta=in1.beta;
    _omega=in1._omega;
    _omega2=in1._omega2;
    _mu=in1._mu;

    _Q=in1._Q;

    _dFdPhase_T1=in1._dFdPhase_T1;
    _dFdPhase_T2=in1._dFdPhase_T2;
    _dFdPhase=in1._dFdPhase;


    _dThetadt_T1=in1._dThetadt_T1;
    _dThetadt_T2=in1._dThetadt_T2;
    _dThetadt=in1._dThetadt;
    return *this;
 }

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass>
void ChenYunACTOriEnergy<FDClass, FDAngleClass>::Calc_Simples(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      
      // _absDeltaTheta_2(j, i) = _Theta->Dx(j,i) * _Theta->Dx(j,i) + _Theta->Dy(j,i) * _Theta->Dy(j,i);
      // _absDeltaTheta(j, i) = sqrt(_absDeltaTheta_2(j, i));
      // if (_absDeltaTheta > _minDelTheta)
      // {
      //   _R(j, i) = 1.0 / _absDeltaTheta(j, i);
      // } else {
      //   _R(j, i) = 1.0 / _minDelTheta;
      // }
      
      // _R3(j, i) = _R(j, i) * _R(j, i) * _R(j, i);

      _Q(j,i) = 1.0 + (1.0 - _mu/_omega) * exp(-_beta * _omega * _Theta->Mag(j,i));

    }
}

template <class FDClass, class FDAngleClass>
void ChenYunACTOriEnergy<FDClass, FDAngleClass>::Calc_dThetadt_T1(){
  double t2 = 0.0;
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      
      _dThetadt_T1(j,i) = _Theta->Mag(j,i) * _Theta->Mag(j,i) * (_Theta->Dxx(j,i) + _Theta->Dyy(j,i));
      _dThetadt_T1(j,i) -= _Theta->Dxy(j,i) * (_Theta->Dx(j,i) + _Theta->Dy(j,i)) +
        _Theta->Dx(j,i) * _Theta->Dxx(j,i) + _Theta->Dy(j,i) * _Theta->Dyy(j,i);
      _dThetadt_T1(j,i) *= (_Phi->F2(j,i)) * _Theta->R3(j,i);

      double tmp = _Theta->Mag(j,i);
      if (abs(tmp) <= 1.0e-12)
      {
        tmp = 1.0e-8;
      }

      t2 = _Phi->Dx(j,i) * _Theta->Dx(j,i) + _Phi->Dy(j,i) * _Theta->Dy(j,i);
      t2 *= 2.0 * (_Phi->F(j,i)) / tmp;

      _dThetadt_T1(j,i) += t2;
      _dThetadt_T1(j,i) *= _alpha;
    }
}


template <class FDClass, class FDAngleClass>
void ChenYunACTOriEnergy<FDClass, FDAngleClass>::Calc_dThetadt_T2(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      
      _dThetadt_T2(j,i) = 2.0 * _Phi->F(j,i) * (
          _Phi->Dx(j,i) * _Theta->Dx(j,i) +
          _Phi->Dy(j,i) * _Theta->Dy(j,i) 
        );

      _dThetadt_T2(j,i) += _Phi->F2(j,i) * (_Theta->Dxx(j,i) + _Theta->Dyy(j,i));
      _dThetadt_T2(j,i) *= _omega2;
    }
}




template <class FDClass, class FDAngleClass>
void ChenYunACTOriEnergy<FDClass, FDAngleClass>::Calc_dThetadt(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      
      _dThetadt(j, i) = _dThetadt_T1(j,i) + _dThetadt_T2(j,i);
      _dThetadt(j,i) /= _Q(j,i); 

    }
}

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass>
void ChenYunACTOriEnergy<FDClass, FDAngleClass>::Calc_dFdPhase(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _dFdPhase(j,i)= -2.0 * _Phi->F(j,i) * _Theta->Mag(j,i);
      _dFdPhase(j,i) -= _omega2 * _Phi->F(j,i) * _Theta->Mag(j,i) * _Theta->Mag(j,i);
    }
}

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass>
void ChenYunACTOriEnergy<FDClass, FDAngleClass>::Calc_All(){
  Calc_Simples();
  Calc_dThetadt_T1();
  Calc_dThetadt_T2();
  Calc_dThetadt();

  Calc_dFdPhase();

}


#endif
