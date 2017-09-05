#ifndef TAKANGLE_H
#define TAKANGLE_H

#include "JFDMpiAngle.h"

template <class FDClass>
class TakAngle {
public:
  TakAngle(JMpi inJMpi, const JMat &in_F, const double &inMinMag);
  TakAngle(JMpi inJMpi, const double &inMinMag);
	TakAngle<FDClass> & operator= (const TakAngle<FDClass> &in1); //Write to operator

  //=========================
  // Useful Functions
  // void SetF(const JMat &in_F); // DANGER DON'T USE

  void Calc_FD();
  void Calc_Rs();
  void Calc_All();

  void Update_Theta(const double &dThetadt, const double &dt, const int &y, const int &x);
  //=========================
  // GetPointers
  JMat * FP(){return &_F;};
  FDClass * DP(){return &_D;};

  //=========================
  // GetPointers
  JMat * DxP(){return _D.DxP();};
  JMat * DyP(){return _D.DyP();};

  JMat * DxxP(){return _D.DxxP();};
  JMat * DyyP(){return _D.DyyP();};

  JMat * DxyP(){return _D.DxyP();};
  JMat * D2P(){return _D.D2P();};

  // GetValues
  double F(const int &y, const int &x){return _F.Value(y,x);};

  double Mag(const int &y, const int &x){return _Mag.Value(y,x);};
  double R(const int &y, const int &x){return _R.Value(y,x);};
  double R3(const int &y, const int &x){return _R3.Value(y,x);};

  double Dx(const int &y, const int &x){return _D.Dx(y,x);};
  double Dy(const int &y, const int &x){return _D.Dy(y,x);};
  double Dxy(const int &y, const int &x){return _D.Dxy(y,x);};
  double Dxx(const int &y, const int &x){return _D.Dxx(y,x);};
  double Dyy(const int &y, const int &x){return _D.Dyy(y,x);};
  double D2(const int &y, const int &x){return _D.D2(y,x);};

protected:
  JMpi _MpiObj;
	int _NY, _NX, _Ny;

  JMat _F;
  FDClass _D;
  JMat _Mag , _R, _R3;
  double _MinMag;

  MPI_Status _status;

private:
  void WrapToPi(double &ValIn);
  const double pidouble=2*M_PI;
};

//**********************************************************************
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
void TakAngle<FDClass>::Update_Theta(const double &dThetadt, const double &dt, const int &y, const int &x){
  _F(y,x)+=dThetadt*dt;
  WrapToPi(_F(y,x));
}

#endif
