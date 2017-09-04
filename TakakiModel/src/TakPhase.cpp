#include "TakPhase.h"

//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
TakPhase<FDClass>::TakPhase(JMpi inJMpi,const JMat &in_F) :
  _MpiObj(inJMpi), _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()),
  _Ny(_MpiObj.NYLo()), _F(in_F), _F2(_NY,_NX), _F3(_NY,_NX), _F4(_NY,_NX),
  _F5(_NY,_NX), _D(_MpiObj, _F), _P(_NY,_NX), _dP(_NY,_NX) {}
