#ifndef CYOri_DThetaDT_H
#define CYOri_DThetaDT_H

#include <ChenYunOriLHS_Q.h>
#include <ChenYunOriRHS_Term1.h>
#include <ChenYunOriRHS_Term2.h>

template <class FDClass, class FDAngleClass>
class CYOri_DThetaDT {

public:
    CYOri_DThetaDT(
        CYOriLHS_Q<FDAngleClass>* Q,
        CYOriRHS_Term1<FDClass, FDAngleClass>* RHS_T1,
        CYOriRHS_Term2<FDClass, FDAngleClass>* RHS_T2,
        const double& tau_theta,
        JMpi inJMpi
    );

    CYOri_DThetaDT<FDClass, FDAngleClass>& operator=(
        const CYOri_DThetaDT<FDClass, FDAngleClass>& obj
    );

    double val(const int &y, const int &x) {return _val(y,x);};
    JMat* val_ptr() {return &(_val);};
protected:
    double _tau_theta;
    double _inv_tau_theta;

    CYOriLHS_Q<FDAngleClass>* _Q;
    CYOriRHS_Term1<FDClass, FDAngleClass>* _RHS_T1;
    CYOriRHS_Term2<FDClass, FDAngleClass>* _RHS_T2;

    JMpi _MpiObj;
    int _NY, _NX, _Ny;

    JMat _val;

protected:
    void calc_all();

};


template <class FDClass, class FDAngleClass>
CYOri_DThetaDT<FDClass, FDAngleClass>::CYOri_DThetaDT(
    CYOriLHS_Q<FDAngleClass>* Q;
    CYOriRHS_Term1<FDClass, FDAngleClass>* RHS_T1;
    CYOriRHS_Term2<FDClass, FDAngleClass>* RHS_T2;
    const double& tau_theta,
    JMpi inJMpi
  ) :  _tau_theta(tau_theta), _inv_tau_theta(1.0 / tau_theta),
  _Q(Q), _RHS_T1(RHS_T1), _RHS_T2(RHS_T2),
  _MpiObj(inJMpi),
  _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()), _Ny(_MpiObj.NYLo()),
  _val(_NY, _NX)
{}


template <class FDClass, class FDAngleClass>
CYOri_DThetaDT<FDClass, FDAngleClass>& CYOri_DThetaDT<FDClass, FDAngleClass>::operator=(
    const CYOri_DThetaDT<FDClass, FDAngleClass>& obj)
{
    _tau_theta = obj._tau_theta;
    _inv_tau_theta = obj._inv_tau_theta;

    _Q = obj._Q;
    _RHS_T1 = obj._RHS_T1;
    _RHS_T2 = obj._RHS_T2;


    _MpiObj = obj._MpiObj;
    _NY = obj._NY;
    _NX = obj._NX;
    _Ny = obj._Ny;

    _val = obj._val;

    return *this;
}


template <class FDClass, class FDAngleClass>
void CYOri_DThetaDT<FDClass, FDAngleClass>::calc_all()
{

    for (int j=0; j<_Ny; j++)
        for (int i=0; i<_NX; i++){
            _val(j, i) = _RHS_T1->val(j,i) + _RHS_T2->val(j,i);
            _val(j, i) *= _Q->inv(j,i) * _inv_tau_theta;
        }
}


#endif