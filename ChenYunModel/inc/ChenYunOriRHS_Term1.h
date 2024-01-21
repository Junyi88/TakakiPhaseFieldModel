// This is the calculations for the first term of equation 3 of the chen yun papaer
// https://wulixb.iphy.ac.cn/pdf-content/10.7498/aps.58.124.pdf

// $\nabla \cdot \left[ \alpha m(\phi) \frac{\nabla \theta}{|\nabla \theta|}  \right]$
// Notes in here https://www.evernote.com/shard/s157/nl/19266417/d6b2eb07-7447-8624-554a-f9639725e310?title=Phase%20Field%20Model%20Calculations

#ifndef CYOriRHS_Term1_H
#define CYOriRHS_Term1_H

#include "TakPhase.h"
#include "TakAngle.h"


template <class FDClass, class FDAngleClass>
class CYOriRHS_Term1 {

public:
    CYOriRHS_Term1(
        TakPhase<FDClass>* inPhi,
        TakAngle<FDAngleClass>* inTheta,
        const double& alpha,
        const double& min_theta,
        JMpi inJMpi
    );

    CYOriRHS_Term1<FDClass, FDAngleClass>& operator=(
        const CYOriRHS_Term1<FDClass, FDAngleClass>& obj
    );

    double val(const int &y, const int &x) {return _val(y,x);};
    JMat* val_ptr() {return &(_val);};

protected:
    double _alpha;
    double _min_theta;

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

// Definition
template <class FDClass, class FDAngleClass>
CYOriRHS_Term1<FDClass, FDAngleClass>::CYOriRHS_Term1(
    TakPhase<FDClass>* inPhi,
    TakAngle<FDAngleClass>* inTheta,
    const double& alpha,
    const double& min_theta,
    JMpi inJMpi
) :  _alpha(alpha), _min_theta(min_theta),
    _MpiObj(inJMpi),
    _NY(_MpiObj.NYGl()),  _NX(_MpiObj.NX()),  _Ny(_MpiObj.NYLo()),
    _term1(_NY,_NX), _term2(_NY, _NX)
{}


template <class FDClass, class FDAngleClass>
CYOriRHS_Term1<FDClass, FDAngleClass>& CYOriRHS_Term1<FDClass, FDAngleClass>::operator=(
    const CYOriRHS_Term1<FDClass, FDAngleClass>& obj)
{
    _alpha = obj._alpha;
    _min_theta = obj._min_theta;

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
void CYOriRHS_Term1<FDClass, FDAngleClass>::calc_term_1()
{
  double tmp = 0.0;
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      
      tmp = _Theta->Mag(j,i);
      if (abs(tmp) > _min_theta)
      {
        _term1(j,i) = tmp * tmp * (_Theta->Dxx(j,i) + _Theta->Dyy(j,i));
      
        _term1(j,i) -= _Theta->Dxy(j,i) * (_Theta->Dx(j,i) + _Theta->Dy(j,i)) +
            _Theta->Dx(j,i) * _Theta->Dxx(j,i) + _Theta->Dy(j,i) * _Theta->Dyy(j,i);
        _term1(j,i) *= (_Phi->F2(j,i)) * _Theta->R3(j,i);
      } else 
      {
        _term1(j,i) = 0.0;
      }

    }
}



template <class FDClass, class FDAngleClass>
void CYOriRHS_Term1<FDClass, FDAngleClass>::calc_term_2()
{
  double tmp = 0.0;
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      

      tmp = _Theta->Mag(j,i);
      if (abs(tmp) > _min_theta)
      {
        _term2(j,i) = _Phi->Dx(j,i) * _Theta->Dx(j,i) + _Phi->Dy(j,i) * _Theta->Dy(j,i);
        _term2(j,i) *= 2.0 * (_Phi->F(j,i)) / tmp;
      } else 
      {
        _term2(j,i) = 0.0;
      }

    }
}


template <class FDClass, class FDAngleClass>
void CYOriRHS_Term1<FDClass, FDAngleClass>::calc_all()
{
    calc_term_1();
    calc_term_2();

    for (int j=0; j<_Ny; j++)
        for (int i=0; i<_NX; i++){
            _val(j, i) = _term1(j, i) + _term2(j, i);
            _val(j, i) *= _alpha;
        }

}

#endif
