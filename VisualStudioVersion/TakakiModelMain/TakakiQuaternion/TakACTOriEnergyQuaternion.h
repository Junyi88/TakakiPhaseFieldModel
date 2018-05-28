#ifndef TAKACTORIENERGYQUATERNION_H
#define TAKACTORIENERGYQUATERNION_H

#include "TakPhase.h"
#include "TakQuaternion.h"

//##========================================================================
template <class FDClass, class FDAngleClass>
class TakACTOriEnergyQuaternion {
public:
	TakACTOriEnergyQuaternion(TakPhase<FDClass> * inPhi, TakQuaternion<FDAngleClass> * inTheta,
		const double &insConst, const double &inMTheta0, const double &inInvPhiMin,
		JMpi inJMpi);
	TakACTOriEnergyQuaternion<FDClass, FDAngleClass> & operator= (
		const TakACTOriEnergyQuaternion<FDClass, FDAngleClass> &in1); //Write to operator

	void Calc_dFdPhase();

	void Calc_MTheta();

	void Calc_DRs();
	void Calc_qDqs();
	void Calc_D2Qs();

	void Calc_dThetadt();

	void Calc_All();

	// Getter Functions

	double dFdPhase(const int &y, const int &x) { return _dFdPhase(y, x); };
	double dThetadt1(const int &y, const int &x) { return _dThetadt1(y, x); };
	double dThetadt2(const int &y, const int &x) { return _dThetadt2(y, x); };
	double dThetadt3(const int &y, const int &x) { return _dThetadt3(y, x); };
	double dThetadt4(const int &y, const int &x) { return _dThetadt4(y, x); };

	JMat * dFdPhasePointer() { return &(_dFdPhase); };
	JMat * dThetadt1Pointer() { return &(_dThetadt1); };
	JMat * dThetadt2Pointer() { return &(_dThetadt2); };
	JMat * dThetadt3Pointer() { return &(_dThetadt3); };
	JMat * dThetadt4Pointer() { return &(_dThetadt4); };

	JMat * MThetaPointer() { return &(_MTheta); };

protected:
	TakPhase<FDClass> * _Phi;
	TakQuaternion<FDAngleClass> * _Theta;
	JMpi _MpiObj;
	int _NY, _NX, _Ny;

	double _sConst, _MTheta0, _Ms, _s2, _shalf, _InvPhiMin, _InvPhi;
	JMat _dFdPhase;
	JMat _MTheta;
	JMat _dThetadt1, _dThetadt2, _dThetadt3, _dThetadt4;

	JMat _DRX, _DRY;
	JMat _qDqX, _qDqY;
	JMat _DqDq, _qDDq;

//	JMat _Term1, _Term2, _Term3;
};

//**********************************************************************
//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass, class FDAngleClass>
TakACTOriEnergyQuaternion<FDClass, FDAngleClass>::TakACTOriEnergyQuaternion(
	TakPhase<FDClass> * inPhi, TakQuaternion<FDAngleClass> * inTheta,
	const double &insConst, const double &inMTheta0, const double &inInvPhiMin,
	JMpi inJMpi) : _Phi(inPhi), _Theta(inTheta), _MpiObj(inJMpi),
	_NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()), _Ny(_MpiObj.NYLo()),
	_sConst(insConst), _MTheta0(inMTheta0), _Ms(_sConst*_MTheta0), 
	_s2(2.0*_sConst), _shalf(0.5*_sConst),
	_InvPhiMin(inInvPhiMin), _InvPhi(inInvPhiMin),
	_dFdPhase(_NY, _NX), _MTheta(_NY, _NX), 
	_dThetadt1(_NY, _NX), _dThetadt2(_NY, _NX),
	_dThetadt3(_NY, _NX), _dThetadt4(_NY, _NX),
	_DRX(_NY, _NX),_DRY(_NY, _NX),
	_qDqX(_NY, _NX),_qDqY(_NY, _NX),
	_DqDq(_NY, _NX),_qDDq(_NY, _NX) {}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass, class FDAngleClass>
TakACTOriEnergyQuaternion<FDClass, FDAngleClass> & TakACTOriEnergyQuaternion<FDClass, FDAngleClass>::operator=
 (const TakACTOriEnergyQuaternion<FDClass, FDAngleClass> &in1) {

	_Phi = in1._Phi;
	_Theta = in1._Theta;
	_MpiObj = in1._MpiObj;
	_NY = in1._NY;
	_NX = in1._NX;
	_Ny = in1._Ny;

	_sConst = in1._sConst;
	_MTheta0 = in1._MTheta0;
	_Ms = in1._Ms;
	_s2 = in1._s2;
	_shalf = in1._shalf;
	_InvPhiMin = in1._InvPhiMin;
	_dFdPhase = in1._dFdPhase;

	_MTheta = in1._MTheta;
	_dThetadt1 = in1._dThetadt1;
	_dThetadt2 = in1._dThetadt2;
	_dThetadt3 = in1._dThetadt3;
	_dThetadt4 = in1._dThetadt4;


	_DRX = in1._DRX;
	_DRY = in1._DRY;
	_qDqX = in1._qDqX;
	_qDqY = in1._qDqY;
	_DqDq = in1._DqDq;
	_qDDq = in1._qDDq;

//	_Term1 = in1._Term1;
//	_Term2 = in1._Term2;
//	_Term3 = in1._Term3;
	return *this;
}

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakACTOriEnergyQuaternion<FDClass, FDAngleClass>::Calc_dFdPhase() {
	for (int j = 0; j<_Ny; j++)
		for (int i = 0; i<_NX; i++) {
			_dFdPhase(j, i) = _s2*(_Phi->F(j, i))*(_Theta->Mag(j, i));
		}
}

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakACTOriEnergyQuaternion<FDClass, FDAngleClass>::Calc_MTheta() {
	for (int j = 0; j<_Ny; j++)
		for (int i = 0; i<_NX; i++) {
			_MTheta(j, i) = _Ms*(1.0 - (_Phi->P(j, i)));
		}
}

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakACTOriEnergyQuaternion<FDClass, FDAngleClass>::Calc_All() {
	Calc_dFdPhase();
	Calc_MTheta();

	Calc_DRs();
	Calc_qDqs();
	Calc_D2Qs();

	Calc_dThetadt();

}

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakACTOriEnergyQuaternion<FDClass, FDAngleClass>::Calc_DRs() {
	for (int j = 0; j<_Ny; j++)
		for (int i = 0; i<_NX; i++) {
			_DRX(j, i) = (_Theta->Dx1(j, i))*(_Theta->Dxx1(j, i));
			_DRX(j, i) += (_Theta->Dy1(j, i))*(_Theta->Dxy1(j, i));
			_DRX(j, i) += (_Theta->Dx2(j, i))*(_Theta->Dxx2(j, i));
			_DRX(j, i) += (_Theta->Dy2(j, i))*(_Theta->Dxy2(j, i));
			_DRX(j, i) += (_Theta->Dx3(j, i))*(_Theta->Dxx3(j, i));
			_DRX(j, i) += (_Theta->Dy3(j, i))*(_Theta->Dxy3(j, i));
			_DRX(j, i) += (_Theta->Dx4(j, i))*(_Theta->Dxx4(j, i));
			_DRX(j, i) += (_Theta->Dy4(j, i))*(_Theta->Dxy4(j, i));
			_DRX(j, i) *= -_Theta->R3(j,i);

			_DRY(j, i) = (_Theta->Dx1(j, i))*(_Theta->Dxy1(j, i));
			_DRY(j, i) += (_Theta->Dy1(j, i))*(_Theta->Dyy1(j, i));
			_DRY(j, i) += (_Theta->Dx2(j, i))*(_Theta->Dxy2(j, i));
			_DRY(j, i) += (_Theta->Dy2(j, i))*(_Theta->Dyy2(j, i));
			_DRY(j, i) += (_Theta->Dx3(j, i))*(_Theta->Dxy3(j, i));
			_DRY(j, i) += (_Theta->Dy3(j, i))*(_Theta->Dyy3(j, i));
			_DRY(j, i) += (_Theta->Dx4(j, i))*(_Theta->Dxy4(j, i));
			_DRY(j, i) += (_Theta->Dy4(j, i))*(_Theta->Dyy4(j, i));
			_DRY(j, i) *= -_Theta->R3(j, i);
		}

}

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakACTOriEnergyQuaternion<FDClass, FDAngleClass>::Calc_qDqs() {
	for (int j = 0; j<_Ny; j++)
		for (int i = 0; i<_NX; i++) {
			_qDqX(j, i) = (_Theta->F1(j, i))*(_Theta->Dx1(j, i));
			_qDqX(j, i) += (_Theta->F2(j, i))*(_Theta->Dx2(j, i));
			_qDqX(j, i) += (_Theta->F3(j, i))*(_Theta->Dx3(j, i));
			_qDqX(j, i) += (_Theta->F4(j, i))*(_Theta->Dx4(j, i));

			_qDqY(j, i) = (_Theta->F1(j, i))*(_Theta->Dy1(j, i));
			_qDqY(j, i) += (_Theta->F2(j, i))*(_Theta->Dy2(j, i));
			_qDqY(j, i) += (_Theta->F3(j, i))*(_Theta->Dy3(j, i));
			_qDqY(j, i) += (_Theta->F4(j, i))*(_Theta->Dy4(j, i));
		}

}

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakACTOriEnergyQuaternion<FDClass, FDAngleClass>::Calc_D2Qs() {
	for (int j = 0; j<_Ny; j++)
		for (int i = 0; i<_NX; i++) {
			_DqDq(j, i) = (_Theta->Dx1(j, i))*(_Theta->Dx1(j, i));
			_DqDq(j, i) += (_Theta->Dy1(j, i))*(_Theta->Dy1(j, i));
			_DqDq(j, i) += (_Theta->Dx2(j, i))*(_Theta->Dx2(j, i));
			_DqDq(j, i) += (_Theta->Dy2(j, i))*(_Theta->Dy2(j, i));
			_DqDq(j, i) += (_Theta->Dx3(j, i))*(_Theta->Dx3(j, i));
			_DqDq(j, i) += (_Theta->Dy3(j, i))*(_Theta->Dy3(j, i));
			_DqDq(j, i) += (_Theta->Dx4(j, i))*(_Theta->Dx4(j, i));
			_DqDq(j, i) += (_Theta->Dy4(j, i))*(_Theta->Dy4(j, i));

			_qDDq(j, i) = (_Theta->F1(j, i))*(
				_Theta->Dxx1(j, i)+ _Theta->Dyy1(j, i));
			_qDDq(j, i) += (_Theta->F2(j, i))*(
				_Theta->Dxx2(j, i) + _Theta->Dyy2(j, i));
			_qDDq(j, i) += (_Theta->F3(j, i))*(
				_Theta->Dxx3(j, i) + _Theta->Dyy3(j, i));
			_qDDq(j, i) += (_Theta->F4(j, i))*(
				_Theta->Dxx4(j, i) + _Theta->Dyy4(j, i));
		}

}

// @@ -- Function ----------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakACTOriEnergyQuaternion<FDClass, FDAngleClass>::Calc_dThetadt() {
	for (int j = 0; j<_Ny; j++)
		for (int i = 0; i<_NX; i++) {
			
			if ((_Phi->F(j, i)) >= _InvPhiMin)
				_InvPhi=1.0 / (_Phi->F(j, i));
			else
				_InvPhi = 1.0 / (_InvPhiMin);
			
			// Q1 ----- 
			_dThetadt1(j, i) = _Theta->Dxx1(j, i)*_Theta->Dxx1(j, i) +
				_Theta->Dyy1(j, i)*_Theta->Dyy1(j, i);
			_dThetadt1(j, i) -= (_Theta->Dx1(j, i) * _qDqX(j, i) +
				_Theta->Dy1(j, i) * _qDqY(j, i));
			_dThetadt1(j, i) -= _Theta->F1(j, i)*
				(_DqDq(j,i)+_qDDq(j,i));
			_dThetadt1(j, i) *= _Theta->R(j, i);

			_dThetadt1(j, i) += _DRX(j, i)*(_Theta->Dx1(j, i) -
				_Theta->F1(j, i)*_qDqX(j, i));
			_dThetadt1(j, i) += _DRY(j, i)*(_Theta->Dy1(j, i) -
				_Theta->F1(j, i)*_qDqY(j, i));

			_dThetadt1(j, i) *= _MTheta(j, i) *_shalf*_Phi->F(j, i);

			_dThetadt1(j, i) += _MTheta(j, i) * _sConst * _Theta->R(j,i) *(
				_Phi->Dx(j,i) * (_Theta->Dx1(j,i)-_Theta->F1(j,i)*_qDqX(j,i)) + 
				_Phi->Dy(j, i) * (_Theta->Dy1(j, i) - _Theta->F1(j, i)*_qDqY(j, i))
				);

			// Q2 ----- 
			_dThetadt2(j, i) = _Theta->Dxx2(j, i)*_Theta->Dxx2(j, i) +
				_Theta->Dxx2(j, i)*_Theta->Dxx2(j, i);
			_dThetadt2(j, i) -= (_Theta->Dx2(j, i) * _qDqX(j, i) +
				_Theta->Dy2(j, i) * _qDqY(j, i));
			_dThetadt2(j, i) -= _Theta->F2(j, i)*
				(_DqDq(j, i) + _qDDq(j, i));
			_dThetadt2(j, i) *= _Theta->R(j, i);

			_dThetadt2(j, i) += _DRX(j, i)*(_Theta->Dx2(j, i) -
				_Theta->F2(j, i)*_qDqX(j, i));
			_dThetadt2(j, i) += _DRY(j, i)*(_Theta->Dy2(j, i) -
				_Theta->F2(j, i)*_qDqY(j, i));

			_dThetadt2(j, i) *= _MTheta(j, i) *_shalf*_Phi->F(j, i);

			_dThetadt2(j, i) += _MTheta(j, i) * _sConst * _Theta->R(j, i) *(
				_Phi->Dx(j, i) * (_Theta->Dx2(j, i) - _Theta->F2(j, i)*_qDqX(j, i)) +
				_Phi->Dy(j, i) * (_Theta->Dy2(j, i) - _Theta->F2(j, i)*_qDqY(j, i))
				);

			// Q3 ----- 
			_dThetadt3(j, i) = _Theta->Dxx3(j, i)*_Theta->Dxx3(j, i) +
				_Theta->Dxx3(j, i)*_Theta->Dxx3(j, i);
			_dThetadt3(j, i) -= (_Theta->Dx3(j, i) * _qDqX(j, i) +
				_Theta->Dy3(j, i) * _qDqY(j, i));
			_dThetadt3(j, i) -= _Theta->F3(j, i)*
				(_DqDq(j, i) + _qDDq(j, i));
			_dThetadt3(j, i) *= _Theta->R(j, i);

			_dThetadt3(j, i) += _DRX(j, i)*(_Theta->Dx3(j, i) -
				_Theta->F3(j, i)*_qDqX(j, i));
			_dThetadt3(j, i) += _DRY(j, i)*(_Theta->Dy3(j, i) -
				_Theta->F3(j, i)*_qDqY(j, i));

			_dThetadt3(j, i) *= _MTheta(j, i) *_shalf*_Phi->F(j, i);

			_dThetadt3(j, i) += _MTheta(j, i) * _sConst * _Theta->R(j, i) *(
				_Phi->Dx(j, i) * (_Theta->Dx3(j, i) - _Theta->F3(j, i)*_qDqX(j, i)) +
				_Phi->Dy(j, i) * (_Theta->Dy3(j, i) - _Theta->F3(j, i)*_qDqY(j, i))
				);

			// Q4 ----- 
			_dThetadt4(j, i) = _Theta->Dxx4(j, i)*_Theta->Dxx4(j, i) +
				_Theta->Dxx4(j, i)*_Theta->Dxx4(j, i);
			_dThetadt4(j, i) -= (_Theta->Dx4(j, i) * _qDqX(j, i) +
				_Theta->Dy4(j, i) * _qDqY(j, i));
			_dThetadt4(j, i) -= _Theta->F4(j, i)*
				(_DqDq(j, i) + _qDDq(j, i));
			_dThetadt4(j, i) *= _Theta->R(j, i);

			_dThetadt4(j, i) += _DRX(j, i)*(_Theta->Dx4(j, i) -
				_Theta->F4(j, i)*_qDqX(j, i));
			_dThetadt4(j, i) += _DRY(j, i)*(_Theta->Dy4(j, i) -
				_Theta->F4(j, i)*_qDqY(j, i));

			_dThetadt4(j, i) *= _MTheta(j, i) *_shalf*_Phi->F(j, i);

			_dThetadt4(j, i) += _MTheta(j, i) * _sConst * _Theta->R(j, i)* (
				_Phi->Dx(j, i) * (_Theta->Dx4(j, i) - _Theta->F4(j, i)*_qDqX(j, i)) +
				_Phi->Dy(j, i) * (_Theta->Dy4(j, i) - _Theta->F4(j, i)*_qDqY(j, i))
				);
		}

}

#endif