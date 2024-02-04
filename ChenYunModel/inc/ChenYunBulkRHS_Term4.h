// This is the calculations for the first term of equation 3 of the chen yun papaer
// https://wulixb.iphy.ac.cn/pdf-content/10.7498/aps.58.124.pdf


#ifndef CYBulkRHS_Term4_H
#define CYBulkRHS_Term4_H

#include "TakPhase.h"
#include "TakAngle.h"


template <class FDClass, class FDAngleClass>
class CYBulkRHS_Term4 {

public:
    CYBulkRHS_Term4(
        TakPhase<FDClass>* inPhi,
        TakAngle<FDAngleClass>* inTheta,
        const double& omega,
        JMpi inJMpi
    );

    CYBulkRHS_Term4<FDClass, FDAngleClass>& operator=(
        const CYBulkRHS_Term4<FDClass, FDAngleClass>& obj
    );

    double val(const int &y, const int &x) {return _val(y,x);};
    JMat* val_ptr() {return &(_val);};

protected:
    double _omega2;

    TakPhase<FDClass>* _Phi;
    TakAngle<FDAngleClass>* _Theta;
    JMpi _MpiObj;
    int _NY, _NX, _Ny;
    JMat _val;

protected:
    void calc_all();

};


template <class FDClass, class FDAngleClass>
CYBulkRHS_Term4<FDClass, FDAngleClass>::CYBulkRHS_Term4(
    TakPhase<FDClass>* inPhi,
    TakAngle<FDAngleClass>* inTheta,
    const double& omega,
    JMpi inJMpi
) :  _omega2(omega * omega),
    _Phi(inPhi), _Theta(inTheta),
    _MpiObj(inJMpi),
    _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()), _Ny(_MpiObj.NYLo()),
    _val(_NY, _NX)
{}


template <class FDClass, class FDAngleClass>
CYBulkRHS_Term4<FDClass, FDAngleClass>& CYBulkRHS_Term4<FDClass, FDAngleClass>::operator=(
    const CYBulkRHS_Term4<FDClass, FDAngleClass>& obj)
{
    _omega2 = obj._omega2;

    _Phi = obj._Phi;
    _Theta = obj._Theta;

    _MpiObj = obj._MpiObj;
    _NY = obj._NY;
    _NX = obj._NX;
    _Ny = obj._Ny;

    _val = obj._val;

    return *this;
}


template <class FDClass, class FDAngleClass>
void CYBulkRHS_Term4<FDClass, FDAngleClass>::calc_all()
{

    for (int j=0; j<_Ny; j++)
        for (int i=0; i<_NX; i++){
            _val(j, i) = -_omega2 * _Phi->F(j, i) * _Theta->Mag(j, i) * _Theta->Mag(j, i);
        }
}

#endif
