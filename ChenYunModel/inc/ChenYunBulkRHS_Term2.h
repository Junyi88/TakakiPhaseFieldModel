// This is the calculations for the first term of equation 3 of the chen yun papaer
// https://wulixb.iphy.ac.cn/pdf-content/10.7498/aps.58.124.pdf


#ifndef CYBulkRHS_Term2_H
#define CYBulkRHS_Term2_H

#include "TakPhase.h"
#include "TakAngle.h"


template <class FDClass>
class CYBulkRHS_Term2 {

public:
    CYBulkRHS_Term2(
        TakPhase<FDClass>* inPhi,
        const double& epsilon,
        const double& w,
        JMpi inJMpi
    );

    CYBulkRHS_Term2<FDClass>& operator=(
        const CYBulkRHS_Term2<FDClass>& obj
    );

    double val(const int &y, const int &x) {return _val(y,x);};
    JMat* val_ptr() {return &(_val);};

protected:
    double _a2;

    TakPhase<FDClass>* _Phi;
    JMpi _MpiObj;
    int _NY, _NX, _Ny;
    JMat _val;

protected:
    void calc_all();

};


template <class FDClass>
CYBulkRHS_Term2<FDClass>::CYBulkRHS_Term2(
    TakPhase<FDClass>* inPhi,
    const double& epsilon,
    const double& w,
    JMpi inJMpi
) :  _a2(9.0 * epsilon * epsilon / (w *w)),
    _Phi(inPhi),
    _MpiObj(inJMpi),
    _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()), _Ny(_MpiObj.NYLo()),
    _val(_NY, _NX)
{}


template <class FDClass>
CYBulkRHS_Term2<FDClass>& CYBulkRHS_Term2<FDClass>::operator=(
    const CYBulkRHS_Term2<FDClass>& obj)
{
    _a2 = obj._a2;

    _Phi = obj._Phi;

    _MpiObj = obj._MpiObj;
    _NY = obj._NY;
    _NX = obj._NX;
    _Ny = obj._Ny;

    _val = obj._val;

    return *this;
}


template <class FDClass>
void CYBulkRHS_Term2<FDClass>::calc_all()
{

    for (int j=0; j<_Ny; j++)
        for (int i=0; i<_NX; i++){
            _val(j, i) = _a2 * (1.0 - _Phi->F(j, i));
        }
}

#endif
