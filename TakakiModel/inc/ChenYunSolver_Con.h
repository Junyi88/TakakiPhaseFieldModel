#ifndef ChenYunSolver_Con_H
#define ChenYunSolver_Con_H

#include "TakPhase.h"
#include "TakAngle.h"
#include "ChenYunACBulkEnergy.h"
#include "ChenYunACTOriEnergy.h"
#include "BasicChemPotential.h"
#include "TakACChemEnergy.h"

//##========================================================================
template <class FDClass, class FDAngleClass>
class ChenYunSolver_Con {
public:
  ChenYunSolver_Con(TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta, BasicChemPotential<FDConClass> * inCon,
    ChenYunACBulkEnergy<FDClass> * inBulkEnergy, ChenYunACTOriEnergy<FDClass, FDAngleClass> * inOriEnergy,
    TakACChemEnergy<FDConClass> * inChemEnergy,
    const double &inMPhiConst, const double &indt, const double& tauPhi, const double& tauTheta,
    JMpi inJMpi);
	ChenYunSolver_Con<FDClass, FDAngleClass> & operator= (
    const ChenYunSolver_Con<FDClass, FDAngleClass> &in1); //Write to operator


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

  void Update_Con();
  void Update_Con(const double &dtimeCustom);

  // Getter Functions
  JMat * dEtadtPointer() {return &(_dEtadt);};

protected:
  TakPhase<FDClass> * _Phi;
  TakAngle<FDAngleClass> * _Theta;
  BasicChemPotential<FDConClass> * _Con;
  ChenYunACBulkEnergy<FDClass> * _BulkEnergy;
  ChenYunACTOriEnergy<FDClass, FDAngleClass> * _OriEnergy;

  TakACChemEnergy<FDConClass> * _ChemEnergy;

  JMpi _MpiObj;
	int _NY, _NX, _Ny;

  double _MPhiConst;
  double _dt;
  double _tauPhi;
  double _tauTheta;
  JMat _dEtadt;
  // JMat _dThetadt;

};

//***************************************************************************
//##========================================================================
// @@ -- Constructor ----------------------------------------------------
template <class FDClass, class FDAngleClass>
ChenYunSolver_Con<FDClass, FDAngleClass>::ChenYunSolver_Con(TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta,
   BasicChemPotential<FDConClass> * inCon,
   ChenYunACBulkEnergy<FDClass> * inBulkEnergy, ChenYunACTOriEnergy<FDClass, FDAngleClass> * inOriEnergy,
   TakACChemEnergy<FDConClass> * inChemEnergy,
   const double &inMPhiConst, const double &indt, const double& tauPhi, const double& tauTheta,
   JMpi inJMpi) : _Phi(inPhi), _Theta(inTheta), _Con(inCon),
   _BulkEnergy(inBulkEnergy), _OriEnergy(inOriEnergy), _ChemEnergy(inChemEnergy),
   _MpiObj(inJMpi),
   _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()), _Ny(_MpiObj.NYLo()),
   _MPhiConst(inMPhiConst), _dt(indt), _tauPhi(tauPhi), _tauTheta(tauTheta),
   _dEtadt(_NY,_NX) {Step_NoUpdate();}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass, class FDAngleClass>
ChenYunSolver_Con<FDClass, FDAngleClass> & ChenYunSolver_Con<FDClass, FDAngleClass>::operator= (
 const ChenYunSolver_Con<FDClass, FDAngleClass> &in1) {

  _Phi=in1._Phi;
  _Theta=in1._Theta;
  _Con=in1._Con;
  _BulkEnergy=in1._BulkEnergy;
  _OriEnergy=in1._OriEnergy;
  _ChemEnergy=in1._ChemEnergy;
  _MpiObj=in1._MpiObj;
  _NY=in1._NY;
  _NX=in1._NX;
  _Ny=in1._Ny;
  _MPhiConst=in1._MPhiConst;
  _dt=in1._dt;
  _tauPhi=in1.tauPhi;
  _tauTheta=in1.tauTheta;

  _dEtadt=in1._dEtadt;
  // _dThetadt=in1._dThetadt;

  return *this;
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void ChenYunSolver_Con<FDClass, FDAngleClass>::Step_NoUpdate(){
  _Phi->Calc_All();
  _Theta->Calc_All();
  _Con->Calc_All();
  _BulkEnergy->Calc_All();
  _OriEnergy->Calc_All();
  _ChemEnergy->Calc_All();
  Calc_dEtadt();
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void ChenYunSolver_Con<FDClass, FDAngleClass>::Calc_dEtadt(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _dEtadt(j,i)=(_BulkEnergy->dFdPhase(j,i));
      _dEtadt(j,i)+=(_OriEnergy->dFdPhase(j,i));
      
      double tmp = 0.0;
      double C3 = _Con->con(j,i);
      C3 = C3 * C3 * C3;
      double C4 = C3 * _Con->con(j,i);
      double C5 = C4 * _Con->con(j,i);

      if (( _Con->con(j,i) > 0.0) && ( _Con->con(j,i) < 1.0))
      {
        tmp = 10.0 * C3 - 15.0 * C4 + 6.0 * C5;
        _dEtadt(j,i)*=tmp;
      } else if (_Con->con(j,i) <= 0.0)
      {
        _dEtadt(j, i) *= 0.0;
      }
      
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void ChenYunSolver_Con<FDClass, FDAngleClass>::Update_Eta(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Phi->Update_Eta(_dEtadt(j,i) / _tauPhi,_dt,j,i);
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void ChenYunSolver_Con<FDClass, FDAngleClass>::Update_Theta(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){

      double tmp = 0.0;
      double C3 = _Con->con(j,i);
      C3 = C3 * C3 * C3;
      double C4 = C3 * _Con->con(j,i);
      double C5 = C4 * _Con->con(j,i);

      if (( _Con->con(j,i) > 0.0) && ( _Con->con(j,i) < 1.0))
      {
        tmp = 10.0 * C3 - 15.0 * C4 + 6.0 * C5;
      } else if (_Con->con(j,i) <= 0.0)
      {
        tmp = 0.0;
      }

      _Theta->Update_Theta(tmp * _OriEnergy->dThetadt(j,i) / _tauTheta,_dt,j,i);
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void ChenYunSolver_Con<FDClass, FDAngleClass>::Update_Eta(const double &dtimeCustom){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Phi->Update_Eta(_dEtadt(j,i),dtimeCustom,j,i);
    }
}


// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass, class FDConClass>
void TakakiSolverAngleCon<FDClass, FDAngleClass, FDConClass>::Update_Con(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Con->Update_Con(_ChemEnergy->dcondt()->Value(j,i),_dt,j,i);
    }
}

template <class FDClass, class FDAngleClass, class FDConClass>
void TakakiSolverAngleCon<FDClass, FDAngleClass, FDConClass>::Update_Con(const double &dtimeCustom){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Con->Update_Con(_ChemEnergy->dcondt()->Value(j,i),dtimeCustom,j,i);
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void ChenYunSolver_Con<FDClass, FDAngleClass>::Update_Theta(const double &dtimeCustom){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){

      double tmp = 0.0;
      double C3 = _Con->con(j,i);
      C3 = C3 * C3 * C3;
      double C4 = C3 * _Con->con(j,i);
      double C5 = C4 * _Con->con(j,i);

      if (( _Con->con(j,i) > 0.0) && ( _Con->con(j,i) < 1.0))
      {
        tmp = 10.0 * C3 - 15.0 * C4 + 6.0 * C5;
      } else if (_Con->con(j,i) <= 0.0)
      {
        tmp = 0.0;
      }

      _Theta->Update_Theta(tmp * _OriEnergy->dThetadt(j,i) / _tauTheta ,dtimeCustom,j,i);
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void ChenYunSolver_Con<FDClass, FDAngleClass>::Step_All(){
  Update_Eta();
  Update_Theta();
  Update_Con();
  Step_NoUpdate();
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void ChenYunSolver_Con<FDClass, FDAngleClass>::Step_All(const double &dtimeCustom){
  Update_Eta(dtimeCustom);
  Update_Theta(dtimeCustom);
  Update_Con(dtimeCustom);
  Step_NoUpdate();
}

#endif
