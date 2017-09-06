#ifndef TAKACTORIENERGY_H
#define TAKACTORIENERGY_H

#include "TakPhase.h"
#include "TakAngle.h"

//##========================================================================
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

  JMat * dFdPhasePointer() {return &(_dFdPhase);};
  JMat * dThetadtPointer() {return &(_dThetadt);};
  JMat * MThetaPointer() {return &(_MTheta);};

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

//**********************************************************************
//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass, class FDAngleClass>
TakACTOriEnergy<FDClass, FDAngleClass>::TakACTOriEnergy(
  TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta,
  const double &insConst, const double &inMTheta0,  const double &inInvPhiMin,
  JMpi inJMpi) : _Phi(inPhi), _Theta(inTheta), _MpiObj(inJMpi),
  _NY(_MpiObj.NYGl()),  _NX(_MpiObj.NX()),  _Ny(_MpiObj.NYLo()),
  _sConst(insConst), _MTheta0(inMTheta0), _Ms(_sConst*_MTheta0), _s2(2.0*_sConst),
  _InvPhiMin(inInvPhiMin),
  _dFdPhase(_NY,_NX), _MTheta(_NY,_NX), _dAngleFront(_NY,_NX),
  _dAngleRear(_NY,_NX), _dThetadt(_NY,_NX) {}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass, class FDAngleClass>
TakACTOriEnergy<FDClass, FDAngleClass> & TakACTOriEnergy<FDClass, FDAngleClass>::operator=
  (const TakACTOriEnergy<FDClass, FDAngleClass> &in1) {

    _Phi=in1._Phi;
    _Theta=in1._Theta;
    _MpiObj=in1._MpiObj;
  	_NY=in1._NY;
    _NX=in1._NX;
    _Ny=in1._Ny;

    _sConst=in1._sConst;
    _MTheta0=in1._MTheta0;
    _Ms=in1._Ms;
    _s2=in1._s2;
    _InvPhiMin=in1._InvPhiMin;
    _dFdPhase=in1._dFdPhase;

    _MTheta=in1._MTheta;
    _dAngleFront=in1._dAngleFront;
    _dAngleRear=in1._dAngleRear;
    _dThetadt=in1._dThetadt;
    return *this;
 }

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakACTOriEnergy<FDClass, FDAngleClass>::Calc_dFdPhase(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _dFdPhase(j,i)=_s2*(_Phi->F(j,i))*(_Theta->Mag(j,i));
    }
}

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakACTOriEnergy<FDClass, FDAngleClass>::Calc_MTheta(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _MTheta(j,i)=_Ms*(1.0-(_Phi->P(j,i)));
    }
}

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakACTOriEnergy<FDClass, FDAngleClass>::Calc_dAngleFront(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _dAngleFront(j,i)=(_Theta->Dy(j,i))*(_Theta->Dy(j,i))*(_Theta->Dxx(j,i));
      _dAngleFront(j,i)-=(_Theta->Dy(j,i))*(_Theta->Dxy(j,i))*(_Theta->Dx(j,i));

      _dAngleFront(j,i)+=(_Theta->Dx(j,i))*(_Theta->Dx(j,i))*(_Theta->Dyy(j,i));
      _dAngleFront(j,i)-=(_Theta->Dx(j,i))*(_Theta->Dxy(j,i))*(_Theta->Dy(j,i));

      _dAngleFront(j,i)*=_Theta->R3(j,i);
    }
}

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakACTOriEnergy<FDClass, FDAngleClass>::Calc_dAngleRear(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _dAngleRear(j,i)=(_Phi->Dx(j,i))*(_Theta->Dx(j,i));
      _dAngleRear(j,i)+=(_Phi->Dy(j,i))*(_Theta->Dy(j,i));

      if ((_Phi->F(j,i))>=_InvPhiMin)
        _dAngleRear(j,i)*=2.0*(_Theta->R(j,i))/(_Phi->F(j,i));
      else
        _dAngleRear(j,i)*=2.0*(_Theta->R(j,i))/(_InvPhiMin);
    }
}

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakACTOriEnergy<FDClass, FDAngleClass>::Calc_dThetadt(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _dThetadt(j,i)=_dAngleFront(j,i)+_dAngleRear(j,i);
      _dThetadt(j,i)*=_MTheta(j,i);
    }
}

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakACTOriEnergy<FDClass, FDAngleClass>::Calc_All(){
  Calc_dFdPhase();
  Calc_MTheta();
  Calc_dAngleFront();
  Calc_dAngleRear();
  Calc_dThetadt();

}


#endif
