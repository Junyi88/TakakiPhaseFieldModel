#include "TakPhase.h"

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
      _P(j,i)=(10.0*_F3(j,i))-(15.0*_F4(j,i))+(6.0*_F5(j,i));
    }
};

// @@ -- Functions  ----------------------------------------------------
template <class FDClass>
void TakPhase<FDClass>::Calc_dP(){
   for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _dP(j,i)=(30.0*_F2(j,i))-(60.0*_F3(j,i))+(30.0*_F4(j,i));
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
}
