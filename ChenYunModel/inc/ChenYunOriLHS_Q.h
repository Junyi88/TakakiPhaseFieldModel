#ifndef CYOriLHS_Q_H
#define CYOriLHS_Q_H

#include "TakPhase.h"
#include "TakAngle.h"


template <class FDAngleClass>
class CYOriLHS_Q {

public:
    CYOriLHS_Q(
        TakAngle<FDAngleClass>* inTheta,
        const double& beta,
        const double& mu,
        const double& omega,
        JMpi inJMpi
    );

    CYOriLHS_Q<FDAngleClass>& operator=(
        const CYOriLHS_Q<FDAngleClass>& obj
    );

    double val(const int &y, const int &x) {return _val(y,x);};
    JMat* val_ptr() {return &(_val);};
    double inv(const int &y, const int &x) {return _inv(y,x);};
    JMat* inv_ptr() {return &(_inv);};

protected:
    double _beta;
    double _mu;
    double _omega;
    double _beta_omega;
    double _mu_omega;

    TakAngle<FDAngleClass>* _Theta;
    JMpi _MpiObj;
    int _NY, _NX, _Ny;

    JMat _exponent;
    JMat _val;
    JMat _inv;

public:
    void calc_exponent();
    void calc_all();

};


template <class FDAngleClass>
CYOriLHS_Q<FDAngleClass>::CYOriLHS_Q(
    TakAngle<FDAngleClass>* inTheta,
    const double& beta,
    const double& mu,
    const double& omega,
    JMpi inJMpi
) :  _beta(beta),
    _mu(mu),
    _omega(omega),
    _beta_omega(beta*omega),
    _mu_omega(mu/omega),
    _Theta(inTheta),
    _MpiObj(inJMpi),
    _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()), _Ny(_MpiObj.NYLo()),
    _exponent(_NY,_NX), _val(_NY, _NX), _inv(_NY, _NX)
{}


template <class FDAngleClass>
CYOriLHS_Q<FDAngleClass>& CYOriLHS_Q<FDAngleClass>::operator=(
    const CYOriLHS_Q<FDAngleClass>& obj)
{
    _beta = obj._beta;
    _mu = obj._mu;
    _omega = obj._omega;
    _beta_omega = obj._beta_omega;
    _mu_omega = obj._mu_omega;

    _Theta = obj._Theta;
    _MpiObj = obj._MpiObj;
    _NY = obj._NY;
    _NX = obj._NX;
    _Ny = obj._Ny;

    _exponent = obj._exponent;
    _val = obj._val;
    _inv = obj._inv;

    return *this;
}


template <class FDAngleClass>
void CYOriLHS_Q<FDAngleClass>::calc_exponent()
{

  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
        _exponent(j,i) = exp(-_beta_omega * _Theta->Mag(j, i));
    }
}


template <class FDAngleClass>
void CYOriLHS_Q<FDAngleClass>::calc_all()
{
  calc_exponent();

  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
        _val(j,i) = 1.0 - _exponent(j,i) + _mu_omega * _exponent(j,i);
        _inv(j,i) = 1.0 / _val(j,i);
    }
}

#endif
