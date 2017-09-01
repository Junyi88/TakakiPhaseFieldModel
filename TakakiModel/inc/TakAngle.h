#ifndef TAKANGLE_H
#define TAKANGLE_H

#include "JFDMpiAngle.h"

template <class FDClass>
class TakAngle {
public:
  TakAngle(JMpi inJMpi, const JMat &in_F, const double &inMinMag);
  TakAngle(JMpi inJMpi, const double &inMinMag);
	TakAngle<FDClass> & operator= (const TakAngle<FDClass> &in1); //Write to operator

  //=========================
  // Useful Functions
  // void SetF(const JMat &in_F); // DANGER DON'T USE

  void Calc_FD();
  void Calc_Rs();
  void Calc_All();

  //=========================
  // GetPointers
  JMat * FP(){return &_F;};
  FDClass * DP(){return &_D;};

  // GetValues
  double F(const int &y, const int &x){return _F.Value(y,x);};

  double Mag(const int &y, const int &x){return _Mag.Value(y,x);};
  double R(const int &y, const int &x){return _R.Value(y,x);};
  double R3(const int &y, const int &x){return _R3.Value(y,x);};

  double Dx(const int &y, const int &x){return _D.Dx(y,x);};
  double Dy(const int &y, const int &x){return _D.Dy(y,x);};
  double Dxy(const int &y, const int &x){return _D.Dxy(y,x);};
  double Dxx(const int &y, const int &x){return _D.Dxx(y,x);};
  double Dyy(const int &y, const int &x){return _D.Dyy(y,x);};
  double D2(const int &y, const int &x){return _D.D2(y,x);};

protected:
  JMpi _MpiObj;
	int _NY, _NX, _Ny;

  JMat _F;
  FDClass _D;
  JMat _Mag , _R, _R3;
  double _MinMag;

  MPI_Status _status;
};


#endif
