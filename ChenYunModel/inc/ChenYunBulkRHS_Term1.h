// This is the calculations for the first term of equation 3 of the chen yun papaer
// https://wulixb.iphy.ac.cn/pdf-content/10.7498/aps.58.124.pdf


#ifndef CYBulkRHS_Term1_H
#define CYBulkRHS_Term1_H

#include "TakPhase.h"
#include "TakAngle.h"


template <class FDClass>
class CYBulkRHS_Term1 {

public:
    CYBulkRHS_Term1(
        TakPhase<FDClass>* inPhi,
        const double& epsilon,
        JMpi inJMpi
    );

    CYBulkRHS_Term1<FDClass>& operator=(
        const CYBulkRHS_Term1<FDClass>& obj
    );

    double val(const int &y, const int &x) {return _val(y,x);};
    JMat* val_ptr() {return &(_val);};

protected:
    double _epsilon2;

    TakPhase<FDClass>* _Phi;
    JMpi _MpiObj;
    int _NY, _NX, _Ny;
    JMat _val;

protected:
    void calc_all();

};


template <class FDClass>
CYBulkRHS_Term1<FDClass>::CYBulkRHS_Term1(
    TakPhase<FDClass>* inPhi,
    const double& epsilon,
    JMpi inJMpi
) :  _epsilon2(epsilon * epsilon),
    _Phi(inPhi),
    _MpiObj(inJMpi),
    _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()), _Ny(_MpiObj.NYLo()),
    _val(_NY, _NX)
{}


template <class FDClass>
CYBulkRHS_Term1<FDClass>& CYBulkRHS_Term1<FDClass>::operator=(
    const CYBulkRHS_Term1<FDClass>& obj)
{
    _epsilon2 = obj._epsilon2;

    _Phi = obj._Phi;

    _MpiObj = obj._MpiObj;
    _NY = obj._NY;
    _NX = obj._NX;
    _Ny = obj._Ny;

    _val = obj._val;

    return *this;
}


template <class FDClass>
void CYBulkRHS_Term1<FDClass>::calc_all()
{

    for (int j=0; j<_Ny; j++)
        for (int i=0; i<_NX; i++){
            _val(j, i) = _epsilon2 * _Phi->D2(j,i);
        }
}

#endif
