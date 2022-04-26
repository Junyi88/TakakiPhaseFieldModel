template <class FDClass>
BasicChemPotential::BasicChemPotential(JMpi inJMpi, const double& kappa, const JMat& con):
  _MpiObj(inJMpi), _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()),
  _Ny(_MpiObj.NYLo()), 
  kappa_(kappa), con_(con), mu_(_NY, _NX),
  D_con_(_MpiObj, con_), D_mu_(_MpiObj, mu_);
{
  Calc_All();
}

template <class FDClass>
BasicChemPotential::BasicChemPotential(JMpi inJMpi, const double& kappa):
  _MpiObj(inJMpi), _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()),
  _Ny(_MpiObj.NYLo()), 
  kappa_(kappa), con_(_NY,_NX), mu_(_NY, _NX),
  D_con_(_MpiObj, con_), D_mu_(_MpiObj, mu_);
{
  Calc_All();
}

template <class FDClass>
void BasicChemPotential::Calc_Diff_Con()
{
  D_con_.TransferAll();
  D_con_.Calc_All();
}

template <class FDClass>
void BasicChemPotential::Calc_Mu()
{
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++)
    {
      mu_.Replace(j, i, -kappa * D_con_.D2(j, i));
    }
}

template <class FDClass>
void BasicChemPotential::Calc_Diff_Mu()
{
  D_mu_.TransferAll();
  D_mu_.Calc_All();  
}

template <class FDClass>
void BasicChemPotential::Calc_All()
{
  Calc_Diff_Con();
  Calc_Mu();
  Calc_Diff_Mu();
}
