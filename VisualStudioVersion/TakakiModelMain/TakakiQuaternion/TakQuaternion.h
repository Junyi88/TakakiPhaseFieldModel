#ifndef TAKQUARTERNION_H
#define TAKQUARTERNION_H


#include "JFDMpi.h"
template <class FDClass>
class TakQuaternion {
public:
	TakQuaternion(JMpi inJMpi, 
		const JMat &in_F1, const JMat &in_F2,
		const JMat &in_F3, const JMat &in_F4,
		const double &inMinMag);
	TakQuaternion(JMpi inJMpi, const double &inMinMag);
	TakQuaternion<FDClass> & operator= (const TakQuaternion<FDClass> &in1); //Write to operator

	//=========================
	// Useful Functions

	void Calc_FD();
	void Calc_Rs();
	void Calc_All();
	void Update_Theta(
		const double &dF1dt, 
		const double &dF2dt,
		const double &dF3dt,
		const double &dF4dt,
		const double &dt, const int &y, const int &x);
	//=========================
	// GetPointers
	JMat * FP1() { return &_F1; };
	FDClass * DP1() { return &_D1; };
	JMat * FP2() { return &_F2; };
	FDClass * DP2() { return &_D2; };
	JMat * FP3() { return &_F3; };
	FDClass * DP3() { return &_D3; };
	JMat * FP4() { return &_F4; };
	FDClass * DP4() { return &_D4; };

	//=========================
	// GetPointers
	JMat * DxP1() { return _D1.DxP(); };
	JMat * DyP1() { return _D1.DyP(); };
	JMat * DxxP1() { return _D1.DxxP(); };
	JMat * DyyP1() { return _D1.DyyP(); };
	JMat * DxyP1() { return _D1.DxyP(); };
	JMat * D2P1() { return _D1.D2P(); };

	JMat * DxP2() { return _D2.DxP(); };
	JMat * DyP2() { return _D2.DyP(); };
	JMat * DxxP2() { return _D2.DxxP(); };
	JMat * DyyP2() { return _D2.DyyP(); };
	JMat * DxyP2() { return _D2.DxyP(); };
	JMat * D2P2() { return _D2.D2P(); };

	JMat * DxP3() { return _D3.DxP(); };
	JMat * DyP3() { return _D3.DyP(); };
	JMat * DxxP3() { return _D3.DxxP(); };
	JMat * DyyP3() { return _D3.DyyP(); };
	JMat * DxyP3() { return _D3.DxyP(); };
	JMat * D2P3() { return _D3.D2P(); };

	JMat * DxP4() { return _D4.DxP(); };
	JMat * DyP4() { return _D4.DyP(); };
	JMat * DxxP4() { return _D4.DxxP(); };
	JMat * DyyP4() { return _D4.DyyP(); };
	JMat * DxyP4() { return _D4.DxyP(); };
	JMat * D2P4() { return _D4.D2P(); };

	// GetValues
	double F(const int &y, const int &x) { return _F.Value(y, x); };

	double Mag(const int &y, const int &x) { return _Mag.Value(y, x); };
	double R(const int &y, const int &x) { return _R.Value(y, x); };
	double R3(const int &y, const int &x) { return _R3.Value(y, x); };

	double Dx1(const int &y, const int &x) { return _D1.Dx(y, x); };
	double Dy1(const int &y, const int &x) { return _D1.Dy(y, x); };
	double Dxy1(const int &y, const int &x) { return _D1.Dxy(y, x); };
	double Dxx1(const int &y, const int &x) { return _D1.Dxx(y, x); };
	double Dyy1(const int &y, const int &x) { return _D1.Dyy(y, x); };
	double D21(const int &y, const int &x) { return _D1.D2(y, x); };

	double Dx2(const int &y, const int &x) { return _D2.Dx(y, x); };
	double Dy2(const int &y, const int &x) { return _D2.Dy(y, x); };
	double Dxy2(const int &y, const int &x) { return _D2.Dxy(y, x); };
	double Dxx2(const int &y, const int &x) { return _D2.Dxx(y, x); };
	double Dyy2(const int &y, const int &x) { return _D2.Dyy(y, x); };
	double D22(const int &y, const int &x) { return _D2.D2(y, x); };

	double Dx3(const int &y, const int &x) { return _D3.Dx(y, x); };
	double Dy3(const int &y, const int &x) { return _D3.Dy(y, x); };
	double Dxy3(const int &y, const int &x) { return _D3.Dxy(y, x); };
	double Dxx3(const int &y, const int &x) { return _D3.Dxx(y, x); };
	double Dyy3(const int &y, const int &x) { return _D3.Dyy(y, x); };
	double D23(const int &y, const int &x) { return _D3.D2(y, x); };

	double Dx4(const int &y, const int &x) { return _D4.Dx(y, x); };
	double Dy4(const int &y, const int &x) { return _D4.Dy(y, x); };
	double Dxy4(const int &y, const int &x) { return _D4.Dxy(y, x); };
	double Dxx4(const int &y, const int &x) { return _D4.Dxx(y, x); };
	double Dyy4(const int &y, const int &x) { return _D4.Dyy(y, x); };
	double D24(const int &y, const int &x) { return _D4.D2(y, x); };

protected:
	JMpi _MpiObj;
	int _NY, _NX, _Ny;

	JMat _F1,_F2,_F3,_F3;
	FDClass _D1,_D2,_D3,_D4;
	JMat _Mag, _R, _R3;
	double _MinMag;

	MPI_Status _status;

private:
	const double pidouble = 2 * M_PI;
};

//**********************************************************************
//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
TakQuaternion<FDClass>::TakQuaternion(JMpi inJMpi, 
	const JMat &in_F1, const JMat &in_F2, 
	const JMat &in_F3, const JMat &in_F4,
	const double &inMinMag) :
	_MpiObj(inJMpi), _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()),
	_Ny(_MpiObj.NYLo()), 
	_F1(in_F1), _F2(in_F2), _F3(in_F3), _F4(in_F4),
	_D1(_MpiObj, _F1), _D2(_MpiObj, _F2), _D2(_MpiObj, _F2), _D4(_MpiObj, _F4),
	_Mag(_NY, _NX), _R(_NY, _NX), _R3(_NY, _NX), _MinMag(inMinMag) {}

// @@ -- Constructor ----------------------------------------------------
template <class FDClass>
TakQuaternion<FDClass>::TakQuaternion(JMpi inJMpi, const double &inMinMag) :
	_MpiObj(inJMpi), _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()),
	_Ny(_MpiObj.NYLo()), 
	_F1(_NY, _NX), _F2(_NY, _NX), _F3(_NY, _NX), _F4(_NY, _NX),
	_D1(_MpiObj, _F1), _D2(_MpiObj, _F2), _D2(_MpiObj, _F2), _D4(_MpiObj, _F4),
	_Mag(_NY, _NX), _R(_NY, _NX), _R3(_NY, _NX), _MinMag(inMinMag) {}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
TakQuaternion<FDClass> & TakQuaternion<FDClass>::operator= (const TakQuaternion<FDClass> &in1) {
	_MpiObj = in1._MpiObj;
	_NY = in1._NY;
	_NX = in1._NX;
	_Ny = in1._Ny;

	_F1 = in1._F1;
	_D1 = _D1.SetF(_F1); //DANGER
	_F2 = in1._F2;
	_D2 = _D2.SetF(_F2); //DANGER
	_F3 = in1._F3;
	_D3 = _D3.SetF(_F3); //DANGER
	_F4 = in1._F4;
	_D4 = _D4.SetF(_F4); //DANGER

	_Mag = in1._Mag;
	_R = in1._R;
	_R3 = in1._R3;

	_MinMag = in1._MinMag;

	return *this;
}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
void TakQuaternion<FDClass>::Calc_FD() {
	_D1.TransferAll();
	_D1.Calc_All();

	_D2.TransferAll();
	_D2.Calc_All();

	_D3.TransferAll();
	_D3.Calc_All();

	_D4.TransferAll();
	_D4.Calc_All();
}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
void TakQuaternion<FDClass>::Calc_All() {
	Calc_FD();
	Calc_Rs();
}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
void TakQuaternion<FDClass>::Update_Theta(
	const double &dF1dt, 
	const double &dF2dt,
	const double &dF3dt,
	const double &dF4dt,
	const double &dt, const int &y, const int &x) {
	_F1(y, x) += dF1dt*dt;
	_F2(y, x) += dF2dt*dt;
	_F3(y, x) += dF3dt*dt;
	_F4(y, x) += dF4dt*dt;
}


// @@ -- Write over operator ----------------------------------------------------
template <class FDClass>
void TakAngle<FDClass>::Calc_Rs() {
	for (int j = 0; j<_Ny; j++)
		for (int i = 0; i<_NX; i++) {
			_Mag(j, i) = sqrt(
				_D1.Dx(j, i)*_D1.Dx(j, i) + _D1.Dy(j, i)*_D1.Dy(j, i)+
				_D2.Dx(j, i)*_D2.Dx(j, i) + _D2.Dy(j, i)*_D2.Dy(j, i) +
				_D3.Dx(j, i)*_D3.Dx(j, i) + _D3.Dy(j, i)*_D3.Dy(j, i) +
				_D4.Dx(j, i)*_D4.Dx(j, i) + _D4.Dy(j, i)*_D4.Dy(j, i) 
			);

			if (_Mag(j, i) >= _MinMag) {
				_R(j, i) = 1.0 / _Mag(j, i);
				_R3(j, i) = _R(j, i)*_R(j, i);
				_R3(j, i) *= _R(j, i);
			}
			else {
				_R(j, i) = 1.0 / _MinMag;
				_R3(j, i) = _R(j, i)*_R(j, i);
				_R3(j, i) *= _R(j, i);
			}
		}

}

#endif


