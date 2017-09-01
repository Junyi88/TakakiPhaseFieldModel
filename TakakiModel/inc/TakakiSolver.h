#ifndef TAKAKISOLVER_H
#define TAKAKISOLVER_H

#include "TakPhase.h"
#include "TakAngle.h"
#include "TakACBulkEnergy.h"
#include "TakACWallEnergy.h"
#include "TakACGradEnergy.h"
#include "TakACTOriEnergy.h"

//##========================================================================
template <class FDClass, class FDAngleClass>
class TakakiSolver {
public:
  TakakiSolver(TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta,
    const double &insConst, const double &inMTheta0,  const double &inInvPhiMin,
    JMpi inJMpi);
	TakakiSolver<FDClass, FDAngleClass> & operator= (
    const TakakiSolver<FDClass, FDAngleClass> &in1); //Write to operator

  void Calc_dEtadt();
  void Calc_dThetadt();
  void Calc_All();


  void Update_Eta();
  void Update_Theta();
  void Update_Eta(const double &dtimeCustom);
  void Update_Theta(const double &dtimeCustom);

  void First_Step();
  void Step_all(const double &dtimeCustom);
  void Step_all();
  
  // Getter Functions

protected:
  TakPhase<FDClass> * _Phi;
  TakAngle<FDAngleClass> * _Theta;
  TakACBulkEnergy<FDClass> * _BulkEnergy;
  TakACWallEnergy<FDClass> * _WallEnergy;
  TakACGradEnergy<FDClass> * _GradEnergy;
  TakACTOriEnergy<FDClass, FDAngleClass> * _OriEnergy;

  JMpi _MpiObj;
	int _NY, _NX, _Ny;

  double _MPhiConst;
  double _dt;
  JMat _dEtadt;
  JMat _dThetadt;

};

#endif
