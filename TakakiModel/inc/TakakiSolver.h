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
  JMat * dEtadtPointer() {return &(_dEtadt);};

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

//***************************************************************************
//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass, class FDAngleClass>
TakakiSolver<FDClass, FDAngleClass>::TakakiSolver(TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta,
   TakACBulkEnergy<FDClass> * inBulkEnergy, TakACWallEnergy<FDClass> * inWallEnergy,
   TakACGradEnergy<FDClass> * inGradEnergy, TakACTOriEnergy<FDClass, FDAngleClass> * inOriEnergy,
   const double &inMPhiConst, const double &indt,
   JMpi inJMpi) : _Phi(inPhi), _Theta(inTheta),
   _BulkEnergy(inBulkEnergy), _WallEnergy(inWallEnergy), _GradEnergy(inGradEnergy),
   _OriEnergy(inOriEnergy), _MpiObj(inJMpi),
   _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()), _Ny(_MpiObj.NYLo()),
   _MPhiConst(inMPhiConst), _dt(indt),
   _dEtadt(_NY,_NX) {Step_NoUpdate();}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass, class FDAngleClass>
TakakiSolver<FDClass, FDAngleClass> & TakakiSolver<FDClass, FDAngleClass>::operator= (
 const TakakiSolver<FDClass, FDAngleClass> &in1) {

  _Phi=in1._Phi;
  _Theta=in1._Theta;
  _BulkEnergy=in1._BulkEnergy;
  _WallEnergy=in1._WallEnergy;
  _GradEnergy=in1._GradEnergy;
  _OriEnergy=in1._OriEnergy;
  _MpiObj=in1._MpiObj;
  _NY=in1._NY;
  _NX=in1._NX;
  _Ny=in1._Ny;
  _MPhiConst=in1._MPhiConst;
  _dt=in1._dt;
  _dEtadt=in1._dEtadt;
  // _dThetadt=in1._dThetadt;
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolver<FDClass, FDAngleClass>::Step_NoUpdate(){
  _BulkEnergy->Calc_All();
  _WallEnergy->Calc_All();
  _GradEnergy->Calc_All();
  _OriEnergy->Calc_All();
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolver<FDClass, FDAngleClass>::Calc_dEtadt(){
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
void TakakiSolver<FDClass, FDAngleClass>::Update_Eta(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Phi.Update_Eta(_dEtadt(j,i),_dt,j,i);
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolver<FDClass, FDAngleClass>::Update_Theta(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Theta.Update_Theta(_OriEnergy->dThetadt(j,i),_dt,j,i);
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolver<FDClass, FDAngleClass>::Update_Eta(const double &dtimeCustom){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Phi.Update_Eta(_dEtadt(j,i),dtimeCustom,j,i);
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolver<FDClass, FDAngleClass>::Update_Theta(const double &dtimeCustom){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Theta.Update_Theta(_OriEnergy->dThetadt(j,i),dtimeCustom,j,i);
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolver<FDClass, FDAngleClass>::Step_All(){
  Update_Eta();
  Update_Theta();
  Step_NoUpdate();
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolver<FDClass, FDAngleClass>::Step_All(const double &dtimeCustom){
  Update_Eta(dtimeCustom);
  Update_Theta(dtimeCustom);
  Step_NoUpdate();
}

#endif
