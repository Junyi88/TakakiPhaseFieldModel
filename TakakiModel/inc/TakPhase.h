#ifndef TAKPHASE_H
#define TAKPHASE_H

#include "JFDMpi.h"

template <class FDClass>
class TakPhase {
public:
  TakPhase(JMpi inJMpi, const JMat &in_F);
  TakPhase(JMpi inJMpi);
	TakPhase<FDClass> & operator= (const TakPhase<FDClass> &in1); //Write to operator

  //=========================
  // Useful Functions
  // void SetF(const JMat &in_F); // DANGER DON'T USE

  void Calc_FD();
  void Calc_Powers();
  void Calc_P();
  void Calc_dP();
  void Calc_All();

  //=========================
  // GetPointers
  JMat * FP(){return &_F;};
  JMat * F2P(){return &_F2;};
  JMat * PP(){return &_P;};
  JMat * dPP(){return &_dP;};
  FDClass * DP(){return &_D;};

  // GetValues
  double F(const int &y, const int &x){return _F.Value(j,i);};
  double F2(const int &y, const int &x){return _F.Value(j,i);};
  double P(const int &y, const int &x){return _P.Value(j,i);};
  double dP(const int &y, const int &x){return _dP.Value(j,i);};

  double Dx(const int &y, const int &x){return _D.Dx(j,i);};
  double Dy(const int &y, const int &x){return _D.Dy(j,i);};
  double Dxy(const int &y, const int &x){return _D.Dxy(j,i);};
  double Dxx(const int &y, const int &x){return _D.Dxx(j,i);};
  double Dyy(const int &y, const int &x){return _D.Dyy(j,i);};
  double D2(const int &y, const int &x){return _D.D2(j,i);};

protected:
  JMpi _MpiObj;
	int _NY, _NX, _Ny;

  JMat _F;
  JMat _F2, _F3, _F4, _F5;
  FDClass _D;
  JMat _P;
  JMat _dP;

  MPI_Status _status;
};


#endif
