#ifndef TAKAKISOLVERQUATERNIONCON_H
#define TAKAKISOLVERQUATERNIONCON_H

#include "TakPhase.h"
#include "TakQuaternion.h"
#include "BasicChemPotential.h"

#include "TakACBulkEnergy.h"
#include "TakACWallEnergy.h"
#include "TakACGradEnergy.h"
#include "TakACTOriEnergyQuaternion.h"

#include "TakACChemEnergy.h"

//##========================================================================
template <class FDClass, class FDAngleClass>
class TakakiSolverQuaternionCon {
public:
	TakakiSolverQuaternionCon(TakPhase<FDClass> * inPhi, TakQuaternion<FDAngleClass> * inTheta,
		BasicChemPotential<FDClass>* con,
		TakACBulkEnergy<FDClass> * inBulkEnergy, TakACWallEnergy<FDClass> * inWallEnergy,
		TakACGradEnergy<FDClass> * inGradEnergy, TakACTOriEnergyQuaternion<FDClass, FDAngleClass> * inOriEnergy,
		TakACChemEnergy<FDClass>* ChemEnergy,
		const double &inMPhiConst, const double &indt,
		JMpi inJMpi);
	TakakiSolverQuaternionCon<FDClass, FDAngleClass> & operator= (
		const TakakiSolverQuaternionCon<FDClass, FDAngleClass> &in1); //Write to operator


	void Step_NoUpdate();

	void Calc_dEtadt();
	// void Calc_dThetadt();
	//void Calc_All();

	void Update_Eta();
	void Update_Theta();
	void Update_Eta(const double &dtimeCustom);
	void Update_Theta(const double &dtimeCustom);

	void Update_Con();
	void Update_Con(const double &dtimeCustom);

	void Step_All(const double &dtimeCustom);
	void Step_All();

	// Getter Functions
	JMat * dEtadtPointer() { return &(_dEtadt); }
	JMat* dcondt();

protected:
	TakPhase<FDClass> * _Phi;
	TakQuaternion<FDAngleClass> * _Theta;
	BasicChemPotential<FDClass>* _con;

	TakACBulkEnergy<FDClass> * _BulkEnergy;
	TakACWallEnergy<FDClass> * _WallEnergy;
	TakACGradEnergy<FDClass> * _GradEnergy;
	TakACTOriEnergyQuaternion<FDClass, FDAngleClass> * _OriEnergy;

	JMpi _MpiObj;
	int _NY, _NX, _Ny;
	TakACChemEnergy<FDClass>* _ChemEnergy;

	double _MPhiConst;
	double _dt;
	JMat _dEtadt;
	// JMat _dThetadt;
	// JMat _dcondt;

};

#include "TakakiSolverQuaternionCon.template.cpp"

#endif

