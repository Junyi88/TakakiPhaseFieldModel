#ifndef CYBulk_DPhiDT_H
#define CYBulk_DPhiDT_H

#include <ChenYunBulkRHS_Term1.h>
#include <ChenYunBulkRHS_Term2.h>
#include <ChenYunBulkRHS_Term3.h>
#include <ChenYunBulkRHS_Term4.h>

template <class FDClass, class FDAngleClass>
class CYBulk_DPhiDT {

public:
    CYBulk_DPhiDT(
        CYBulkRHS_Term1<FDClass>* RHS_T1,
        CYBulkRHS_Term2<FDClass>* RHS_T2,
        CYBulkRHS_Term3<FDClass, FDAngleClass>* RHS_T3,
        CYBulkRHS_Term4<FDClass, FDAngleClass>* RHS_T4,
        const double& tau_phi,
        JMpi inJMpi
    );

    CYBulk_DPhiDT<FDClass, FDAngleClass>& operator=(
        const CYBulk_DPhiDT<FDClass, FDAngleClass>& obj
    );

    double val(const int &y, const int &x) {return _val(y,x);};
    JMat* val_ptr() {return &(_val);};
protected:
    double _tau_phi;
    double _inv_tau_phi;


    CYBulkRHS_Term1<FDClass>* _RHS_T1;
    CYBulkRHS_Term2<FDClass>* _RHS_T2;
    CYBulkRHS_Term3<FDClass, FDAngleClass>* _RHS_T3;
    CYBulkRHS_Term4<FDClass, FDAngleClass>* _RHS_T4;

    JMpi _MpiObj;
    int _NY, _NX, _Ny;

    JMat _val;

protected:
    void calc_all();

};


template <class FDClass, class FDAngleClass>
CYBulk_DPhiDT<FDClass, FDAngleClass>::CYBulk_DPhiDT(
    CYBulkRHS_Term1<FDClass>* RHS_T1,
    CYBulkRHS_Term2<FDClass>* RHS_T2,
    CYBulkRHS_Term3<FDClass, FDAngleClass>* RHS_T3,
    CYBulkRHS_Term4<FDClass, FDAngleClass>* RHS_T4,
    const double& tau_phi,
    JMpi inJMpi
  ) :  _tau_phi(tau_phi), _inv_tau_phi(1.0 / tau_phi),
  _RHS_T1(RHS_T1), _RHS_T2(RHS_T2), _RHS_T3(RHS_T3), _RHS_T4(RHS_T4),
  _MpiObj(inJMpi),
  _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()), _Ny(_MpiObj.NYLo()),
  _val(_NY, _NX)
{}


template <class FDClass, class FDAngleClass>
CYBulk_DPhiDT<FDClass, FDAngleClass>& CYBulk_DPhiDT<FDClass, FDAngleClass>::operator=(
    const CYBulk_DPhiDT<FDClass, FDAngleClass>& obj)
{
    _tau_phi = obj._tau_phi;
    _inv_tau_phi = obj._inv_tau_phi;

    _RHS_T1 = obj._RHS_T1;
    _RHS_T2 = obj._RHS_T2;
    _RHS_T3 = obj._RHS_T3;
    _RHS_T4 = obj._RHS_T4;

    _MpiObj = obj._MpiObj;
    _NY = obj._NY;
    _NX = obj._NX;
    _Ny = obj._Ny;

    _val = obj._val;

    return *this;
}


template <class FDClass, class FDAngleClass>
void CYBulk_DPhiDT<FDClass, FDAngleClass>::calc_all()
{

    for (int j=0; j<_Ny; j++)
        for (int i=0; i<_NX; i++){
            _val(j, i) = _RHS_T1->val(j,i) + _RHS_T2->val(j,i) + _RHS_T3->val(j,i) + _RHS_T4->val(j,i);
            _val(j, i) *= _inv_tau_phi;
        }
}


#endif