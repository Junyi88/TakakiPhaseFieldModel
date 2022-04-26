#ifndef BASIC_CHEMPOTENTIAL_H
#define BASIC_CHEMPOTENTIAL_H

#include <JFDMpi.h>

//== ChemicalPotential
template <class FDClass>
class BasicChemPotential
{
public:
  BasicChemPotential() = delete;
  BasicChemPotential(const BasicChemPotential& obj) = delete;
  BasicChemPotential& operator= (const BasicChemPotential& obj) = delete;
  ~BasicChemPotential() = default;

  BasicChemPotential(JMpi inJMpi, const double& kappa, const JMat& con);
  BasicChemPotential(JMpi inJMpi, const double& kappa);

  //
  void Calc_Diff_Con();
  void Calc_Mu();
  void Calc_Diff_Mu();
  void Calc_All();
  
  void Update_Con(const double &dcondt, const double &dt, const int &y, const int &x);

  double con(const int &y, const int &x){return con_.Value(y,x);}
  double D2_mu(const int &y, const int &x){return D_mu_(x);};
protected:
  JMpi _MpiObj;
  int _NY, _NX, _Ny;
  double kappa_;

  JMat con_;
  JMat mu_;

  FDClass D_con_;
  FDClass D_mu_;

  MPI_Status _status;
};


#include <BasicChemPotential.template.h>
#endif