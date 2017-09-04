#include "TakPhase.h"

//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
TakPhase<FDClass>::TakPhase(JMpi inJMpi,const JMat &in_F) :
  _MpiObj(inJMpi), _F(_NY,_NX) {}
