#ifndef TAKAKISOLVERANGLECON_H
#define TAKAKISOLVERANGLECON_H

#include "TakPhase.h"
#include "TakAngle.h"
#include "TakACBulkEnergy.h"
#include "TakACWallEnergy.h"
#include "TakACGradEnergy.h"
#include "TakACTOriEnergy.h"

#include "BasicChemPotential.h"
#include "TakACChemEnergy.h"

//##========================================================================
template <class FDClass, class FDAngleClass>
class TakakiSolverAngleCon {
public:
  TakakiSolverAngleCon(TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta, BasicChemPotential<FDClass> * inCon,
    TakACBulkEnergy<FDClass> * inBulkEnergy, TakACWallEnergy<FDClass> * inWallEnergy,
    TakACGradEnergy<FDClass> * inGradEnergy, TakACTOriEnergy<FDClass, FDAngleClass> * inOriEnergy,
    TakACChemEnergy<FDClass> * inChemEnergy,
    const double &inMPhiConst, const double &indt,
    JMpi inJMpi);
	TakakiSolverAngleCon<FDClass, FDAngleClass> & operator= (
    const TakakiSolver<FDClass, FDAngleClass> &in1); //Write to operator


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
  JMat * dEtadtPointer() {return &(_dEtadt);};

protected:
  TakPhase<FDClass> * _Phi;
  TakAngle<FDAngleClass> * _Theta;
  BasicChemPotential<FDClass> * _Con;
  TakACBulkEnergy<FDClass> * _BulkEnergy;
  TakACWallEnergy<FDClass> * _WallEnergy;
  TakACGradEnergy<FDClass> * _GradEnergy;
  TakACTOriEnergy<FDClass, FDAngleClass> * _OriEnergy;
  TakACChemEnergy<FDClass> * _ChemEnergy;

  JMpi _MpiObj;
	int _NY, _NX, _Ny;

  double _MPhiConst;
  double _dt;
  JMat _dEtadt;
  // JMat _dThetadt;

};

//***************************************************************************
//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass, class FDAngleClass>
TakakiSolverAngleCon<FDClass, FDAngleClass>::TakakiSolverAngleCon(TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta,
   BasicChemPotential<FDClass> * inCon,
   TakACBulkEnergy<FDClass> * inBulkEnergy, TakACWallEnergy<FDClass> * inWallEnergy,
   TakACGradEnergy<FDClass> * inGradEnergy, TakACTOriEnergy<FDClass, FDAngleClass> * inOriEnergy,
   TakACChemEnergy<FDClass> * inChemEnergy,
   const double &inMPhiConst, const double &indt,
   JMpi inJMpi) : _Phi(inPhi), _Theta(inTheta), _Con(inCon),
   _BulkEnergy(inBulkEnergy), _WallEnergy(inWallEnergy), _GradEnergy(inGradEnergy),
   _OriEnergy(inOriEnergy), _ChemEnergy(inChemEnergy),
   _MpiObj(inJMpi),
   _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()), _Ny(_MpiObj.NYLo()),
   _MPhiConst(inMPhiConst), _dt(indt),
   _dEtadt(_NY,_NX) {Step_NoUpdate();}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass, class FDAngleClass>
TakakiSolverAngleCon<FDClass, FDAngleClass> & TakakiSolverAngleCon<FDClass, FDAngleClass>::operator= (
 const TakakiSolverAngleCon<FDClass, FDAngleClass> &in1) {

  _Phi=in1._Phi;
  _Theta=in1._Theta;
  _Con=in1._Con;
  _BulkEnergy=in1._BulkEnergy;
  _WallEnergy=in1._WallEnergy;
  _GradEnergy=in1._GradEnergy;
  _OriEnergy=in1._OriEnergy;
  _ChemEnergy=in1._ChemEnergy;
  _MpiObj=in1._MpiObj;
  _NY=in1._NY;
  _NX=in1._NX;
  _Ny=in1._Ny;
  _MPhiConst=in1._MPhiConst;
  _dt=in1._dt;
  _dEtadt=in1._dEtadt;
  // _dThetadt=in1._dThetadt;

  return *this;
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverAngleCon<FDClass, FDAngleClass>::Step_NoUpdate(){
  _Phi->Calc_All();
  _Theta->Calc_All();
  _Con->Calc_All();
  _BulkEnergy->Calc_All();
  _WallEnergy->Calc_All();
  _GradEnergy->Calc_All();
  _OriEnergy->Calc_All();
  _ChemEnergy->Calc_All();
  Calc_dEtadt();
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverAngleCon<FDClass, FDAngleClass>::Calc_dEtadt(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _dEtadt(j,i)=-(_BulkEnergy->dFdPhase(j,i));
      _dEtadt(j,i)-=(_WallEnergy->dFdPhase(j,i));
      _dEtadt(j,i)-=(_GradEnergy->dFdPhase(j,i));
      _dEtadt(j,i)-=(_OriEnergy->dFdPhase(j,i));
      _dEtadt(j,i)*=_MPhiConst;
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverAngleCon<FDClass, FDAngleClass>::Update_Eta(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Phi->Update_Eta(_dEtadt(j,i),_dt,j,i);
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverAngleCon<FDClass, FDAngleClass>::Update_Theta(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Theta->Update_Theta(_OriEnergy->dThetadt(j,i),_dt,j,i);
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverAngleCon<FDClass, FDAngleClass>::Update_Eta(const double &dtimeCustom){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Phi->Update_Eta(_dEtadt(j,i),dtimeCustom,j,i);
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverAngleCon<FDClass, FDAngleClass>::Update_Theta(const double &dtimeCustom){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Theta->Update_Theta(_OriEnergy->dThetadt(j,i),dtimeCustom,j,i);
    }
}


// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverAngleCon<FDClass, FDAngleClass>::Update_Con(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Con->Update_Con(_ChemEnergy->dcondt()->operator(j,i),_dt,j,i);
    }
}

template <class FDClass, class FDAngleClass>
void TakakiSolverAngleCon<FDClass, FDAngleClass>::Update_Con(const double &dtimeCustom){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Con->Update_Con(_ChemEnergy->dcondt()->operator(j,i),dtimeCustom,j,i);
    }
}
// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverAngleCon<FDClass, FDAngleClass>::Step_All(){
  Update_Eta();
  Update_Theta();
  Update_Con();
  Step_NoUpdate();
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverAngleCon<FDClass, FDAngleClass>::Step_All(const double &dtimeCustom){
  Update_Eta(dtimeCustom);
  Update_Theta(dtimeCustom);
  Update_Con(dtimeCustom);
  Step_NoUpdate();
}

#endif
