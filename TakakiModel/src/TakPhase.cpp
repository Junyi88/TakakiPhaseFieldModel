#include "TakPhase.h"

//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
ACBulk<FDClass>::ACBulk(JMpi inJMpi, const JMat &in_F) :
  _MpiObj(inJMpi), _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()),
  _Ny(_MpiObj.NYLo()), _F(in_F), _F2(_NY,_NX), _F3(_NY,_NX), _F4(_NY,_NX),
  _D(_MpiObj, _F) {}

// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
ACBulk<FDClass>::ACBulk(JMpi inJMpi) :
  _MpiObj(inJMpi), _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()),
  _Ny(_MpiObj.NYLo()), _F(_NY,_NX), _F2(_NY,_NX), _F3(_NY,_NX), _F4(_NY,_NX),
  _D(_MpiObj, _F) {}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
ACBulk & ACBulk<FDClass>::operator= (const ACBulk &in1) {
  _MpiObj=in1._MpiObj;
	_NY=in1._NY;
  _NX=in1._NX;
  _Ny=in1._Ny;

  _F=in1._F;
  _F2=in1._F2;
  _F3=in1._F3;
  _F4=in1._F4;
  _D=in1._D;
}
