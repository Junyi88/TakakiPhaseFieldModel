#ifndef TAKPHASE_H
#define TAKPHASE_H

#include "JFDMpi.h"

//## ===================================
template <class FDClass>
class TakPhase {
public:
  TakPhase(JMpi inJMpi,const JMat &in_F);
  TakPhase(JMpi inJMpi);
	TakPhase<FDClass> & operator= (const TakPhase<FDClass> &in1); //Write to operator

  //=========================
  // Useful Functions
  // void SetF(const JMat &in_F); // DANGER DON'T USE

  void Calc_FD();
  void Calc_Powers();
  void Calc_P();
  void Calc_dP();
  void Calc_All();

  void Update_Eta(const double &dEtadt, const double &dt, const int &y, const int &x);
  //=========================
  // GetPointers
  JMat * FP(){return &_F;};
  JMat * F2P(){return &_F2;};
  JMat * F3P(){return &_F3;};
  JMat * F4P(){return &_F4;};
  JMat * F5P(){return &_F5;};

  JMat * PP(){return &_P;};
  JMat * dPP(){return &_dP;};
  FDClass * DP(){return &_D;};

  JMat * DxP(){return _D.DxP();};
  JMat * DyP(){return _D.DyP();};

  JMat * DxxP(){return _D.DxxP();};
  JMat * DyyP(){return _D.DyyP();};

  JMat * DxyP(){return _D.DxyP();};
  JMat * D2P(){return _D.D2P();};

  // GetValues
  double F(const int &y, const int &x){return _F.Value(y,x);};
  double F2(const int &y, const int &x){return _F2.Value(y,x);};
  double F3(const int &y, const int &x){return _F3.Value(y,x);};
  double F4(const int &y, const int &x){return _F4.Value(y,x);};
  double F5(const int &y, const int &x){return _F5.Value(y,x);};

  double P(const int &y, const int &x){return _P.Value(y,x);};
  double dP(const int &y, const int &x){return _dP.Value(y,x);};

  double Dx(const int &y, const int &x){return _D.Dx(y,x);};
  double Dy(const int &y, const int &x){return _D.Dy(y,x);};
  double Dxy(const int &y, const int &x){return _D.Dxy(y,x);};
  double Dxx(const int &y, const int &x){return _D.Dxx(y,x);};
  double Dyy(const int &y, const int &x){return _D.Dyy(y,x);};
  double D2(const int &y, const int &x){return _D.D2(y,x);};

  //


protected:
  JMpi _MpiObj;
	int _NY, _NX, _Ny;

  JMat _F;
  JMat _F2, _F3, _F4, _F5;
  FDClass _D;
  JMat _P;
  JMat _dP;

  MPI_Status _status;
};

//**********************************************************************
//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
TakPhase<FDClass>::TakPhase(JMpi inJMpi,const JMat &in_F) :
  _MpiObj(inJMpi), _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()),
  _Ny(_MpiObj.NYLo()), _F(in_F), _F2(_NY,_NX), _F3(_NY,_NX), _F4(_NY,_NX),
  _F5(_NY,_NX), _D(_MpiObj, _F), _P(_NY,_NX), _dP(_NY,_NX) {}

// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
TakPhase<FDClass>::TakPhase(JMpi inJMpi) :
  _MpiObj(inJMpi), _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()),
  _Ny(_MpiObj.NYLo()), _F(_NY,_NX), _F2(_NY,_NX), _F3(_NY,_NX), _F4(_NY,_NX),
  _F5(_NY,_NX), _D(_MpiObj, _F), _P(_NY,_NX), _dP(_NY,_NX) {}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
TakPhase<FDClass> & TakPhase<FDClass>::operator= (const TakPhase<FDClass> &in1) {
  _MpiObj=in1._MpiObj;
	_NY=in1._NY;
  _NX=in1._NX;
  _Ny=in1._Ny;

  _F=in1._F;
  _F2=in1._F2;
  _F3=in1._F3;
  _F4=in1._F4;
  _F5=in1._F5;

  _D=_D.SetF(_F); //DANGER
  //_D=in1._D;
  _P=in1._P;
  _dP=in1._dP;

  return *this;
}

// @@ -- Functions  ----------------------------------------------------
template <class FDClass>
void TakPhase<FDClass>::Calc_FD(){
   _D.TransferAll();
   _D.Calc_All();
};

// @@ -- Functions  ----------------------------------------------------
template <class FDClass>
void TakPhase<FDClass>::Calc_Powers(){
   for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _F2.Replace(j,i, _F.Value(j,i)*_F.Value(j,i));
      _F3.Replace(j,i, _F2.Value(j,i)*_F.Value(j,i));
      _F4.Replace(j,i, _F3.Value(j,i)*_F.Value(j,i));
      _F5.Replace(j,i ,_F4.Value(j,i)*_F.Value(j,i));
    }

};

// @@ -- Functions  ----------------------------------------------------
template <class FDClass>
void TakPhase<FDClass>::Calc_P(){
   for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){

      if ((_F(j,i)>0.0)&&(_F(j,i)<1.0))
        _P(j,i)=(10.0*_F3(j,i))-(15.0*_F4(j,i))+(6.0*_F5(j,i));
      else if (_F(j,i)<=0.0)
        _P(j,i)=0.0;
      else
        _P(j,i)=1.0;
    }
};

// @@ -- Functions  ----------------------------------------------------
template <class FDClass>
void TakPhase<FDClass>::Calc_dP(){
   for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      if ((_F(j,i)>0.0)&&(_F(j,i)<1.0))
        _dP(j,i)=(30.0*_F2(j,i))-(60.0*_F3(j,i))+(30.0*_F4(j,i));
      else
        _dP(j,i)=0.0;
    }
};

// @@ -- Functions  ----------------------------------------------------
template <class FDClass>
void TakPhase<FDClass>::Calc_All(){
  Calc_Powers();
  Calc_FD();
  Calc_P();
  Calc_dP();
}

// @@ -- Functions  ----------------------------------------------------
template <class FDClass>
void TakPhase<FDClass>::Update_Eta(const double &dEtadt, const double &dt, const int &y, const int &x){
  _F(y,x)+=dEtadt*dt;

  if (_F(y,x) > 1.0) {_F(y,x)  = 1.0;}
  if (_F(y,x) < 0.0) {_F(y,x)  = 0.0;}
}


// ************************************************DKFJDASKFJDKASFJ
#endif
