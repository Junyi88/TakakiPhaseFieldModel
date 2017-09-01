#include "TakPhase.h"

//##========================================================================
// @@ -- Constructor ----------------------------------------------------
ACBulk::ACBulk(JMpi inJMpi, const JMat &in_F) :
  _MpiObj(inJMpi), _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()),
  _Ny(_MpiObj.NYLo()), _F(in_F), _F2(_NY,_NX), _F3(_NY,_NX), _F4(_NY,_NX),
  _D(_MpiObj, in_F) {}
