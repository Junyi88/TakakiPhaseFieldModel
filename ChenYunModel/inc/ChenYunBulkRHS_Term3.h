// This is the calculations for the first term of equation 3 of the chen yun papaer
// https://wulixb.iphy.ac.cn/pdf-content/10.7498/aps.58.124.pdf


#ifndef CYBulkRHS_Term3_H
#define CYBulkRHS_Term3_H

#include "TakPhase.h"
#include "TakAngle.h"


template <class FDClass, class FDAngleClass>
class CYBulkRHS_Term3 {

public:
    CYBulkRHS_Term3(
        TakPhase<FDClass>* inPhi,
        TakAngle<FDAngleClass>* inTheta,
        const double& alpha,
        JMpi inJMpi
    );

    CYBulkRHS_Term3<FDClass, FDAngleClass>& operator=(
        const CYBulkRHS_Term3<FDClass, FDAngleClass>& obj
    );

    double val(const int &y, const int &x) {return _val(y,x);};
    JMat* val_ptr() {return &(_val);};

protected:
    double _alpha_times_2;

    TakPhase<FDClass>* _Phi;
    TakAngle<FDAngleClass>* _Theta;
    JMpi _MpiObj;
    int _NY, _NX, _Ny;
    JMat _val;

public:
    void calc_all();

};


template <class FDClass, class FDAngleClass>
CYBulkRHS_Term3<FDClass, FDAngleClass>::CYBulkRHS_Term3(
    TakPhase<FDClass>* inPhi,
    TakAngle<FDAngleClass>* inTheta,
    const double& alpha,
    JMpi inJMpi
) :  _alpha_times_2(2.0 * alpha),
    _Phi(inPhi), _Theta(inTheta),
    _MpiObj(inJMpi),
    _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()), _Ny(_MpiObj.NYLo()),
    _val(_NY, _NX)
{}


template <class FDClass, class FDAngleClass>
CYBulkRHS_Term3<FDClass, FDAngleClass>& CYBulkRHS_Term3<FDClass, FDAngleClass>::operator=(
    const CYBulkRHS_Term3<FDClass, FDAngleClass>& obj)
{
    _alpha_times_2 = obj._alpha_times_2;

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
void CYBulkRHS_Term3<FDClass, FDAngleClass>::calc_all()
{

    for (int j=0; j<_Ny; j++)
        for (int i=0; i<_NX; i++){
            _val(j, i) = -_alpha_times_2 * _Phi->F(j, i) * _Theta->Mag(j, i);
        }
}

#endif
