#ifndef TAKACTORIENERGYQUATERNION_H
#define TAKACTORIENERGYQUATERNION_H

#include "TakPhase.h"
#include "TakQuaternion.h"

//##========================================================================
template <class FDClass, class FDAngleClass>
class TakACTOriEnergyQuaternion {
public:
	TakACTOriEnergyQuaternion(TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta,
		const double &insConst, const double &inMTheta0, const double &inInvPhiMin,
		JMpi inJMpi);
	TakACTOriEnergyQuaternion<FDClass, FDAngleClass> & operator= (
		const TakACTOriEnergyQuaternion<FDClass, FDAngleClass> &in1); //Write to operator

	void Calc_dFdPhase();

	void Calc_MTheta();
	void Calc_dAngleFront();
	void Calc_dAngleRear();
	void Calc_dThetadt();

	void Calc_All();

	// Getter Functions

	double dFdPhase(const int &y, const int &x) { return _dFdPhase(y, x); };
	double dThetadt(const int &y, const int &x) { return _dThetadt(y, x); };

	JMat * dFdPhasePointer() { return &(_dFdPhase); };
	JMat * dThetadtPointer() { return &(_dThetadt); };
	JMat * MThetaPointer() { return &(_MTheta); };

protected:
	TakPhase<FDClass> * _Phi;
	TakAngle<FDAngleClass> * _Theta;
	JMpi _MpiObj;
	int _NY, _NX, _Ny;

	double _sConst, _MTheta0, _Ms, _s2, _InvPhiMin;
	JMat _dFdPhase;

	JMat _MTheta;
	JMat _dAngleFront;
	JMat _dAngleRear;
	JMat _dThetadt;
};



#endif