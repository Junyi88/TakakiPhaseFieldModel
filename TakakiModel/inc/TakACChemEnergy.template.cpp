template <class FDClass>
TakACChemEnergy<FDClass>::TakACChemEnergy(BasicChemPotential<FDClass>* con,
  const double& MChem, JMpi inJMpi): _con(con), _MpiObj(inJMpi),
  _NY(_MpiObj.NYGl()),  _NX(_MpiObj.NX()), _Ny(_MpiObj.NYLo()),
  _MChem(MChem), _dcondt(_NY, _NX)
{
  Calc_All();
}

template <class FDClass>
TakACChemEnergy<FDClass>& TakACChemEnergy<FDClass>::operator= (const TakACChemEnergy<FDClass>& obj)
{
  _con = obj._con;
  _MpiObj = obj._MpiObj;
  _NY = obj._NY;
  _NX = obj._NX;
  _Ny = obj._Ny;

  _MChem = obj._MChem;
  _dcondt = obj._dcondt;

  return *this;
}

template <class FDClass>
void TakACChemEnergy<FDClass>::Calc_All()
{
  for (unsigned int j = 0; j < _Ny; j++)
  {
    for (unsigned int i = 0; i < _NX; i++)
    {
      _dcondt(j, i) = _MChem * _con->D2_mu(j, i);
    }
  }
}
