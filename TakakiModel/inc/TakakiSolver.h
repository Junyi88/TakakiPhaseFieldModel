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
    TakACBulkEnergy<FDClass> * inBulkEnergy, TakACWallEnergy<FDClass> * inWallEnergy,
    TakACGradEnergy<FDClass> * inGradEnergy, TakACTOriEnergy<FDClass, FDAngleClass> * inOriEnergy,
    const double &inMPhiConst, const double &indt,
    JMpi inJMpi);
	TakakiSolver<FDClass, FDAngleClass> & operator= (
    const TakakiSolver<FDClass, FDAngleClass> &in1); //Write to operator


  void Step_NoUpdate();

  void Calc_dEtadt();
  // void Calc_dThetadt();
  //void Calc_All();

  void Update_Eta();
  void Update_Theta();
  void Update_Eta(const double &dtimeCustom);
  void Update_Theta(const double &dtimeCustom);


  void Step_All(const double &dtimeCustom);
  void Step_All();

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
  // JMat _dThetadt;

};

#endif
