#ifndef TAKACCHEMENERGY_H
#define TAKACCHEMENERGY_H

#include <JFDMpi.h>
#include <JMAT.h>
#include "BasicChemPotential.h"

//== 
template <class FDClass>
class TakACChemEnergy
{
public:
  TakACChemEnergy(BasicChemPotential<FDClass>* con,
    const double& MChem, JMpi inJMpi);

  TakACChemEnergy<FDClass>& operator= (const TakACChemEnergy<FDClass>& obj);
  JMat* dcondt() {return &(_dcondt);}

  void Calc_All();

private:
  BasicChemPotential<FDClass>* _con;

  JMPi _MpiObj;
  int _NY, _NX, _Ny;

  const double _MChem;
  JMat _dcondt;
};

#include <TakACChemEnergy.template.cpp>

#endif