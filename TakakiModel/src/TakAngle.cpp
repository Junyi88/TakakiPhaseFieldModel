#include "TakAngle.h"

//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
TakAngle<FDClass>::TakAngle(JMpi inJMpi, const JMat &in_F, const double &inMinMag) :
  _MpiObj(inJMpi), _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()),
  _Ny(_MpiObj.NYLo()), _F(in_F), _D(_MpiObj, _F),
  _Mag(_NY,_NX), _R(_NY,_NX), _R3(_NY,_NX), _MinMag(inMinMag) {}

// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
TakAngle<FDClass>::TakAngle(JMpi inJMpi, const double &inMinMag) :
  _MpiObj(inJMpi), _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()),
  _Ny(_MpiObj.NYLo()), _F(_NY,_NX), _D(_MpiObj, _F),
  _Mag(_NY,_NX), _R(_NY,_NX), _R3(_NY,_NX), _MinMag(inMinMag) {}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
TakAngle<FDClass> & TakAngle<FDClass>::operator= (const TakAngle<FDClass> &in1) {
  _MpiObj=in1._MpiObj;
	_NY=in1._NY;
  _NX=in1._NX;
  _Ny=in1._Ny;

  _F=in1._F;
  _D=_D.SetF(_F); //DANGER

  _Mag=in1._Mag;
  _R=in1._R;
  _R3=in1._R3;

  _MinMag=in1._MinMag;
}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
void TakAngle<FDClass>::Calc_FD(){
  _D.TransferAll();
  _D.Calc_All();
}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
void TakAngle<FDClass>::Calc_Rs(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Mag(j,i)=sqrt(_D.Dx(j,i)*_D.Dx(j,i)+_D.Dy(j,i)*_D.Dy(j,i));

      if (_Mag(j,i)>=_MinMag){
        _R(j,i)=1.0/_Mag(j,i);
        _R3(j,i)=_R(j,i)*_R(j,i);
        _R3(j,i)*=_R(j,i);
      } else {
        _R(j,i)=1.0/_MinMag;
        _R3(j,i)=_R(j,i)*_R(j,i);
        _R3(j,i)*=_R(j,i);
      }
    }

}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
void TakAngle<FDClass>::Calc_All() {
  Calc_FD();
  Calc_Rs();
}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
void TakAngle<FDClass>::WrapToPi(double &ValIn){
  if (ValIn>M_PI)
    ValIn-=pidouble;
  else if (ValIn<-M_PI)
    ValIn+=pidouble;
}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
double TakAngle<FDClass>::Update_Theta(const double &dThetadt, const double &dt, const int &y, const int &x){
  _F(y,x)+=dThetadt*dt;
  WrapToPi(_F(y,x));
}
