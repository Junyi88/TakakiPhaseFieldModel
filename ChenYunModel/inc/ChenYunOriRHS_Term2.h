// This is the calculations for the first term of equation 3 of the chen yun papaer
// https://wulixb.iphy.ac.cn/pdf-content/10.7498/aps.58.124.pdf

// $\nabla \cdot \left[ \omega^2 n \nabla \theta  \right]$
// Notes in here https://www.evernote.com/shard/s157/nl/19266417/d6b2eb07-7447-8624-554a-f9639725e310?title=Phase%20Field%20Model%20Calculations

#ifndef CYOriRHS_Term2_H
#define CYOriRHS_Term2_H

#include "TakPhase.h"
#include "TakAngle.h"

template <class FDClass, class FDAngleClass>
class CYOriRHS_Term2 {

public:
    CYOriRHS_Term2(
        TakPhase<FDClass>* inPhi,
        TakAngle<FDAngleClass>* inTheta,
        const double& omega,
        JMpi inJMpi
    );

    CYOriRHS_Term2<FDClass, FDAngleClass>& operator=(
        const CYOriRHS_Term2<FDClass, FDAngleClass>& obj
    );

    double val(const int &y, const int &x) {return _val(y,x);};
    JMat* val_ptr() {return &(_val);};

protected:
    double _omega2;

    TakPhase<FDClass>* _Phi;
    TakAngle<FDAngleClass>* _Theta;
    JMpi _MpiObj;
	int _NY, _NX, _Ny;

    JMat _term1;
    JMat _term2;
    JMat _val;

protected:
    void calc_term_1();
    void calc_term_2();
    void calc_all();

};


template <class FDClass, class FDAngleClass>
CYOriRHS_Term2<FDClass, FDAngleClass>::CYOriRHS_Term2(
    TakPhase<FDClass>* inPhi,
    TakAngle<FDAngleClass>* inTheta,
    const double& omega,
    JMpi inJMpi
) :  _omega2(omega * omega),
    _MpiObj(inJMpi),
    _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()), _Ny(_MpiObj.NYLo()),
    _term1(_NY,_NX), _term2(_NY, _NX), _val(_NY, _NX)
{}


template <class FDClass, class FDAngleClass>
CYOriRHS_Term2<FDClass, FDAngleClass>& CYOriRHS_Term2<FDClass, FDAngleClass>::operator=(
    const CYOriRHS_Term2<FDClass, FDAngleClass>& obj)
{
    _omega2 = obj._omega2;

    _Phi = obj._Phi;
    _Theta = obj._Theta;
    _MpiObj = obj._MpiObj;
	_NY = obj._NY;
    _NX = obj._NX;
    _Ny = obj._Ny;

    _term1 = obj._term1;
    _term2 = obj._term2;
    _val = obj._val;

    return *this;
}


template <class FDClass, class FDAngleClass>
void CYOriRHS_Term2<FDClass, FDAngleClass>::calc_term_1()
{

  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
        _term1(j,i) = 2.0 * _Phi->F(j,i) * (
            _Phi->Dx(j,i) * _Theta->Dx(j,i) +
            _Phi->Dy(j,i) * _Theta->Dy(j,i) 
        );

    }
}

template <class FDClass, class FDAngleClass>
void CYOriRHS_Term2<FDClass, FDAngleClass>::calc_term_2()
{

  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
        _term2(j,i) = _Phi->F2(j,i) * (_Theta->Dxx(j,i) + _Theta->Dyy(j,i));
    }
}


template <class FDClass, class FDAngleClass>
void CYOriRHS_Term2<FDClass, FDAngleClass>::calc_all()
{
    calc_term_1();
    calc_term_2();

    for (int j=0; j<_Ny; j++)
        for (int i=0; i<_NX; i++){
            _val(j, i) = _term1(j, i) + _term2(j, i);
            _val(j, i) *= _omega2;
        }
}

#endif
