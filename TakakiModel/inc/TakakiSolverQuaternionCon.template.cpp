template <class FDClass, class FDAngleClass>
TakakiSolverQuaternionCon<FDClass, FDAngleClass>::TakakiSolverQuaternionCon(
  TakPhase<FDClass> * inPhi, TakQuaternion<FDAngleClass> * inTheta, BasicChemPotential<FDClass>* con,
  TakACBulkEnergy<FDClass> * inBulkEnergy, TakACWallEnergy<FDClass> * inWallEnergy,
  TakACGradEnergy<FDClass> * inGradEnergy,
  TakACTOriEnergyQuaternion<FDClass, FDAngleClass> * inOriEnergy,
  TakACChemEnergy<FDClass>* ChemEnergy,
  const double &inMPhiConst, const double &indt,
  JMpi inJMpi) : _Phi(inPhi), _Theta(inTheta), _con(con),
  _BulkEnergy(inBulkEnergy), _WallEnergy(inWallEnergy), _GradEnergy(inGradEnergy),
  _OriEnergy(inOriEnergy), _MpiObj(inJMpi),
  _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()), _Ny(_MpiObj.NYLo()),
  _ChemEnergy(ChemEnergy),
  _MPhiConst(inMPhiConst), _dt(indt),
  _dEtadt(_NY, _NX)
{
  Step_NoUpdate();
}

// @@ -- Write over operator ----------------------------------------------------
template <class FDClass, class FDAngleClass>
TakakiSolverQuaternionCon<FDClass, FDAngleClass>& TakakiSolverQuaternionCon<FDClass, FDAngleClass>::operator= (
  const TakakiSolverQuaternionCon<FDClass, FDAngleClass> &in1) {

  _Phi = in1._Phi;
  _Theta = in1._Theta;
  _con = in1._con;

  _BulkEnergy = in1._BulkEnergy;
  _WallEnergy = in1._WallEnergy;
  _GradEnergy = in1._GradEnergy;
  _OriEnergy = in1._OriEnergy;
  _ChemEnergy = in1._ChemEnergy;

  _MpiObj = in1._MpiObj;
  _NY = in1._NY;
  _NX = in1._NX;
  _Ny = in1._Ny;
  _MPhiConst = in1._MPhiConst;
  _dt = in1._dt;
  _dEtadt = in1._dEtadt;
  // _dThetadt=in1._dThetadt;

  return *this;
}

//-- Step_NoUpdate ----
template <class FDClass, class FDAngleClass>
void TakakiSolverQuaternionCon<FDClass, FDAngleClass>::Step_NoUpdate() {
  _Phi->Calc_All();
  _Theta->Calc_All();
  _con->Calc_All();

  _BulkEnergy->Calc_All();
  _WallEnergy->Calc_All();
  _GradEnergy->Calc_All();
  _OriEnergy->Calc_All();
  Calc_dEtadt();

  _ChemEnergy->Calc_All();
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverQuaternionCon<FDClass, FDAngleClass>::Calc_dEtadt() {
  for (int j = 0; j<_Ny; j++)
    for (int i = 0; i<_NX; i++) {
      _dEtadt(j, i) = -(_BulkEnergy->dFdPhase(j, i));
      _dEtadt(j, i) -= (_WallEnergy->dFdPhase(j, i));
      _dEtadt(j, i) -= (_GradEnergy->dFdPhase(j, i));
      _dEtadt(j, i) -= (_OriEnergy->dFdPhase(j, i));
      _dEtadt(j, i) *= _MPhiConst;
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverQuaternionCon<FDClass, FDAngleClass>::Update_Eta() {
  for (int j = 0; j<_Ny; j++)
    for (int i = 0; i<_NX; i++) {
      _Phi->Update_Eta(_dEtadt(j, i), _dt, j, i);
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverQuaternionCon<FDClass, FDAngleClass>::Update_Theta() {
  for (int j = 0; j<_Ny; j++)
    for (int i = 0; i<_NX; i++) {
      _Theta->Update_Theta(_OriEnergy->dThetadt1(j, i), 
        _OriEnergy->dThetadt2(j, i),
        _OriEnergy->dThetadt3(j, i),
        _OriEnergy->dThetadt4(j, i),
        _dt, j, i);
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverQuaternionCon<FDClass, FDAngleClass>::Update_Eta(const double &dtimeCustom) {
  for (int j = 0; j<_Ny; j++)
    for (int i = 0; i<_NX; i++) {
      _Phi->Update_Eta(_dEtadt(j, i), dtimeCustom, j, i);
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverQuaternionCon<FDClass, FDAngleClass>::Update_Theta(const double &dtimeCustom) {
  for (int j = 0; j<_Ny; j++)
    for (int i = 0; i<_NX; i++) {
      _Theta->Update_Theta(_OriEnergy->dThetadt(j, i), dtimeCustom, j, i);
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverQuaternionCon<FDClass, FDAngleClass>::Update_Con() {
  for (int j = 0; j<_Ny; j++)
    for (int i = 0; i<_NX; i++) {
      _con->Update_Con(_ChemEnergy->dcondt()->Value(j, i), _dt, j, i);
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverQuaternionCon<FDClass, FDAngleClass>::Update_Con(const double &dtimeCustom) {
  for (int j = 0; j<_Ny; j++)
    for (int i = 0; i<_NX; i++) {
      _con->Update_Con(_ChemEnergy->dcondt()->Value(j, i), dtimeCustom, j, i);
    }
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverQuaternionCon<FDClass, FDAngleClass>::Step_All() {
  Update_Eta();
  Update_Theta();
  Update_Con();
  Step_NoUpdate();
}

// @@ ------------------------------------------------------
template <class FDClass, class FDAngleClass>
void TakakiSolverQuaternionCon<FDClass, FDAngleClass>::Step_All(const double &dtimeCustom) {
  Update_Eta(dtimeCustom);
  Update_Theta(dtimeCustom);
  Update_Con(dtimeCustom);
  Step_NoUpdate();
}

