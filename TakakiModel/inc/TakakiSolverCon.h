#ifndef TAKAKISOLVERCON_H
#define TAKAKISOLVERCON_H

#include "TakPhase.h"
#include "TakAngle.h"
#include "BasicChemPotential.h"

#include "TakACBulkEnergy.h"
#include "TakACWallEnergy.h"
#include "TakACGradEnergy.h"
#include "TakACTOriEnergy.h"
#include "TakACChemEnergy.h"
//##========================================================================
template <class FDClass, class FDAngleClass>
class TakakiSolverCon {
public:
  TakakiSolverCon(TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta,
    TakACBulkEnergy<FDClass> * inBulkEnergy, TakACWallEnergy<FDClass> * inWallEnergy,
    TakACGradEnergy<FDClass> * inGradEnergy, TakACTOriEnergy<FDClass, FDAngleClass> * inOriEnergy,
    const double &inMPhiConst, const double &indt,
    JMpi inJMpi);
	TakakiSolverCon<FDClass, FDAngleClass> & operator= (
    const TakakiSolverCon<FDClass, FDAngleClass> &in1); //Write to operator


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
  JMat * dEtadtPointer() {return &(_dEtadt);};

protected:
  TakPhase<FDClass> * _Phi;
  TakAngle<FDAngleClass> * _Theta;
  BasicChemPotential<FDClass>* _con;

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
