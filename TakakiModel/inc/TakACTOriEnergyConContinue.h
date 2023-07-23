#ifndef TAKACTORIENERGYCON_H
#define TAKACTORIENERGYCON_H

#include "TakPhase.h"
#include "TakAngle.h"
#include "BasicChemPotential.h"

//##========================================================================
template <class FDClass, class FDAngleClass, class FDConClass>
class TakACTOriEnergyCon {
public:
  TakACTOriEnergyCon(TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta,
    BasicChemPotential<FDConClass> * inCon,
    const double &insConst, const double &inMTheta0,  const double &inInvPhiMin,
    JMpi inJMpi);
  TakACTOriEnergyCon(TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta,
    BasicChemPotential<FDConClass> * inCon,
    const double &insConst, const double &inMTheta0,  const double &inInvPhiMin,
    JMpi inJMpi, double MThetaMin);
	TakACTOriEnergyCon<FDClass, FDAngleClass, FDConClass> & operator= (const TakACTOriEnergyCon<FDClass, FDAngleClass, FDConClass> &in1); //Write to operator

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
  JMat * dAngleFrontPointer() {return &(_dAngleFront);};
  JMat * dAngleRearPointer() {return &(_dAngleRear);};

protected:
  TakPhase<FDClass> * _Phi;
  TakAngle<FDAngleClass> * _Theta;
  BasicChemPotential<FDConClass> * _Con;
  JMpi _MpiObj;
	int _NY, _NX, _Ny;

  double _sConst, _MTheta0, _Ms, _s2, _InvPhiMin, _MThetaMin;
  JMat _dFdPhase;

  JMat _MTheta;
  JMat _dAngleFront;
  JMat _dAngleRear;
  JMat _dThetadt;
};

//**********************************************************************
//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass, class FDAngleClass, class FDConClass>
TakACTOriEnergyCon<FDClass, FDAngleClass, FDConClass>::TakACTOriEnergyCon(
  TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta, BasicChemPotential<FDConClass> * inCon,
  const double &insConst, const double &inMTheta0,  const double &inInvPhiMin,
  JMpi inJMpi) : _Phi(inPhi), _Theta(inTheta), _Con(inCon),
  _MpiObj(inJMpi),
  _NY(_MpiObj.NYGl()),  _NX(_MpiObj.NX()),  _Ny(_MpiObj.NYLo()),
  _sConst(insConst), _MTheta0(inMTheta0), _Ms(_sConst*_MTheta0), _s2(2.0*_sConst),
  _InvPhiMin(inInvPhiMin), _MThetaMin(0.0),
  _dFdPhase(_NY,_NX), _MTheta(_NY,_NX), _dAngleFront(_NY,_NX),
  _dAngleRear(_NY,_NX), _dThetadt(_NY,_NX) {}

template <class FDClass, class FDAngleClass, class FDConClass>
TakACTOriEnergyCon<FDClass, FDAngleClass, FDConClass>::TakACTOriEnergyCon(
  TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta, BasicChemPotential<FDConClass> * inCon,
  const double &insConst, const double &inMTheta0,  const double &inInvPhiMin,
  JMpi inJMpi, double MThetaMin) : _Phi(inPhi), _Theta(inTheta), _Con(inCon),
  _MpiObj(inJMpi),
  _NY(_MpiObj.NYGl()),  _NX(_MpiObj.NX()),  _Ny(_MpiObj.NYLo()),
  _sConst(insConst), _MTheta0(inMTheta0), _Ms(_sConst*_MTheta0), _s2(2.0*_sConst),
  _InvPhiMin(inInvPhiMin), _MThetaMin(MThetaMin),
  _dFdPhase(_NY,_NX), _MTheta(_NY,_NX), _dAngleFront(_NY,_NX),
  _dAngleRear(_NY,_NX), _dThetadt(_NY,_NX) {}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass, class FDAngleClass, class FDConClass>
TakACTOriEnergyCon<FDClass, FDAngleClass, FDConClass> & TakACTOriEnergyCon<FDClass, FDAngleClass, FDConClass>::operator=
  (const TakACTOriEnergyCon<FDClass, FDAngleClass, FDConClass> &in1) {

    _Phi=in1._Phi;
    _Theta=in1._Theta;
    _Con=in1._Con;
    _MpiObj=in1._MpiObj;
  	_NY=in1._NY;
    _NX=in1._NX;
    _Ny=in1._Ny;

    _sConst=in1._sConst;
    _MTheta0=in1._MTheta0;
    _Ms=in1._Ms;
    _s2=in1._s2;
    _InvPhiMin=in1._InvPhiMin;
    _MThetaMin = in1._MThetaMin;
    _dFdPhase=in1._dFdPhase;

    _MTheta=in1._MTheta;
    _dAngleFront=in1._dAngleFront;
    _dAngleRear=in1._dAngleRear;
    _dThetadt=in1._dThetadt;
    return *this;
 }

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass, class FDConClass>
void TakACTOriEnergyCon<FDClass, FDAngleClass, FDConClass>::Calc_dFdPhase(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _dFdPhase(j,i)=_s2*(_Phi->F(j,i))*(_Theta->Mag(j,i));
    }
}

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass, class FDConClass>
void TakACTOriEnergyCon<FDClass, FDAngleClass, FDConClass>::Calc_MTheta(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){

      if (_Phi->P(j,i) <= 0.98)
      {
        _MTheta(j,i)=_Ms*(1.0-(_Phi->P(j,i)));
      } else {

        double grad = (_MThetaMin - _Ms*0.02) / 0.02;
        _MTheta(j,i)= grad * (_Phi->P(j,i) - 0.98) + _Ms*0.02;

      }
      



    }
}

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass, class FDConClass>
void TakACTOriEnergyCon<FDClass, FDAngleClass, FDConClass>::Calc_dAngleFront(){
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
template <class FDClass, class FDAngleClass, class FDConClass>
void TakACTOriEnergyCon<FDClass, FDAngleClass, FDConClass>::Calc_dAngleRear(){
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
template <class FDClass, class FDAngleClass, class FDConClass>
void TakACTOriEnergyCon<FDClass, FDAngleClass, FDConClass>::Calc_dThetadt(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _dThetadt(j,i)=_dAngleFront(j,i)+_dAngleRear(j,i);
      _dThetadt(j,i)*=_MTheta(j,i);


      // double cc = _Con->con(j,i);
      // if (cc <= 0.0)
      // {
      //   _dThetadt(j,i) = 0.0;
      // } else
      // if (cc < 1.0)
      // {
      //   _dThetadt(j,i) *= cc;
      // }

      double tmp = 0.0;
      double C3 = _Con->con(j,i);
      C3 = C3 * C3 * C3;
      double C4 = C3 * _Con->con(j,i);
      double C5 = C4 * _Con->con(j,i);

      if (( _Con->con(j,i) > 0.0) && ( _Con->con(j,i) < 1.0))
      {
        tmp = 10.0 * C3 - 15.0 * C4 + 6.0 * C5;
        _dThetadt(j,i)*=tmp;
      } else if (_Con->con(j,i) <= 0.0)
      {
        _dThetadt(j, i) = 0.0;
      }


    }
}

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass, class FDConClass>
void TakACTOriEnergyCon<FDClass, FDAngleClass, FDConClass>::Calc_All(){
  Calc_dFdPhase();
  Calc_MTheta();
  Calc_dAngleFront();
  Calc_dAngleRear();
  Calc_dThetadt();

}


#endif
