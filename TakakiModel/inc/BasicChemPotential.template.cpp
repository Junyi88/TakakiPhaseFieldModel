template <class FDClass>
BasicChemPotential<FDClass>::BasicChemPotential(JMpi inJMpi, const double& kappa, const JMat& con):
  _MpiObj(inJMpi), _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()),
  _Ny(_MpiObj.NYLo()), 
  kappa_(kappa), con_(con), mu_(_NY, _NX),
  D_con_(_MpiObj, con_), D_mu_(_MpiObj, mu_)
{
  Calc_All();
}

template <class FDClass>
BasicChemPotential<FDClass>::BasicChemPotential(JMpi inJMpi, const double& kappa):
  _MpiObj(inJMpi), _NY(_MpiObj.NYGl()), _NX(_MpiObj.NX()),
  _Ny(_MpiObj.NYLo()), 
  kappa_(kappa), con_(_NY,_NX), mu_(_NY, _NX),
  D_con_(_MpiObj, con_), D_mu_(_MpiObj, mu_)
{
  Calc_All();
}

template <class FDClass>
void BasicChemPotential<FDClass>::Calc_Diff_Con()
{
  std::cout << "B1" << std::endl;
  D_con_.TransferAll();
  std::cout << "B2" << std::endl;
  D_con_.Calc_All();
  std::cout << "B3" << std::endl;
}

template <class FDClass>
void BasicChemPotential<FDClass>::Calc_Mu()
{
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++)
    {
      mu_.Replace(j, i, -kappa_ * D_con_.D2(j, i));
    }
}

template <class FDClass>
void BasicChemPotential<FDClass>::Calc_Diff_Mu()
{
  D_mu_.TransferAll();
  D_mu_.Calc_All();  
}

template <class FDClass>
void BasicChemPotential<FDClass>::Calc_All()
{
  std::cout << "A1" << std::endl;
  Calc_Diff_Con();
  std::cout << "A2" << std::endl;
  Calc_Mu();
  std::cout << "A3" << std::endl;
  Calc_Diff_Mu();
  std::cout << "A4" << std::endl;
}

template <class FDClass>
void BasicChemPotential<FDClass>::Update_Con(const double &dcondt, const double &dt, const int &y, const int &x)
{

  con_(y, x) += dcondt * dt ;

}
