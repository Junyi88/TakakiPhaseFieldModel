#ifndef TAKPHASE_H
#define TAKPHASE_H

#include "JFDMpi.h"

//## ===================================
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

  void Update_Eta(const double &dEtadt, const double &dt, const int &y, const int &x);
  //=========================
  // GetPointers
  JMat * FP(){return &_F;};
  JMat * F2P(){return &_F2;};
  JMat * PP(){return &_P;};
  JMat * dPP(){return &_dP;};
  FDClass * DP(){return &_D;};

  // GetValues
  double F(const int &y, const int &x){return _F.Value(y,x);};
  double F2(const int &y, const int &x){return _F2.Value(y,x);};
  double F3(const int &y, const int &x){return _F3.Value(y,x);};
  double F4(const int &y, const int &x){return _F4.Value(y,x);};
  double F5(const int &y, const int &x){return _F5.Value(y,x);};

  double P(const int &y, const int &x){return _P.Value(y,x);};
  double dP(const int &y, const int &x){return _dP.Value(y,x);};

  double Dx(const int &y, const int &x){return _D.Dx(y,x);};
  double Dy(const int &y, const int &x){return _D.Dy(y,x);};
  double Dxy(const int &y, const int &x){return _D.Dxy(y,x);};
  double Dxx(const int &y, const int &x){return _D.Dxx(y,x);};
  double Dyy(const int &y, const int &x){return _D.Dyy(y,x);};
  double D2(const int &y, const int &x){return _D.D2(y,x);};

  //


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
