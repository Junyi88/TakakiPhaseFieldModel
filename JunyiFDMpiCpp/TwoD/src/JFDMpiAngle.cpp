#include "JFDMpiAngle.h"

// ## =====================================================================
// @@ -----------------------------------------------
JFDMpi2DReflectHighAngle::JFDMpi2DReflectHighAngle(JMpi inJMpi, JMat &in_F) :
	JFDMpi2DReflectHigh(inJMpi, in_F) {}

// @@ -----------------------------------------------
JFDMpi2DReflectHighAngle & JFDMpi2DReflectHighAngle::operator= (const JFDMpi2DReflectHighAngle &in1){

  _MpiObj=in1._MpiObj;
  _F=in1._F;

  _NY=in1._NY;
  _NX=in1._NX;
  _Ny=in1._Ny;
  _Dx=in1._Dx;
  _Dxx=in1._Dxx;
  _Dy=in1._Dy;
  _Dyy=in1._Dyy;
  _D2=in1._D2;
  _FTop=in1._FTop;
  _FBot=in1._FBot;
  _FTop2=in1._FTop2;
  _FBot2=in1._FBot2;

  _dx=in1._dx;
  _dy=in1._dy;
  _dxdouble=in1._dxdouble;
  _dydouble=in1._dydouble;
  _dxdx=in1._dxdx;
  _dydy=in1._dydy;

  _Dxy=in1._Dxy;
  _FTop3=in1._FTop3;
  _FBot3=in1._FBot3;
  _FTop4=in1._FTop4;
  _FBot4=in1._FBot4;

	return *this;
}

// @@ -----------------------------------------------
double JFDMpi2DReflectHighAngle::WrapAngle(const double &Next, const double &Ref){

  AngleBuffer=Next-Ref;
  if (fabs(AngleBuffer)<=M_PI)
    return AngleBuffer;
  else if ((AngleBuffer)>M_PI)
    return (AngleBuffer-Pi2);
  else if ((AngleBuffer)<-M_PI)
    return (AngleBuffer+Pi2);

  return AngleBuffer;

}

// @@ -----------------------------------------------
void JFDMpi2DReflectHighAngle::CleanUpAllAngles(){

  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _F(j,i)=WrapAngle(_F(j,i),0.0);
    }

}

// @@ ----------------------------------------------------------------
void JFDMpi2DReflectHighAngle::Calc_Dx(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Dx(j,i)=O1C[0]*WrapAngle(FVal(j,i-4),_F(j,i))+
        O1C[1]*WrapAngle(FVal(j,i-3),_F(j,i))+
        O1C[2]*WrapAngle(FVal(j,i-2),_F(j,i))+
        O1C[3]*WrapAngle(FVal(j,i-1),_F(j,i))+
        O1C[5]*WrapAngle(FVal(j,i+1),_F(j,i))+
        O1C[6]*WrapAngle(FVal(j,i+2),_F(j,i))+
        O1C[7]*WrapAngle(FVal(j,i+3),_F(j,i))+
        O1C[8]*WrapAngle(FVal(j,i+4),_F(j,i));
      _Dx(j,i)/=_dx;
    }
}

// @@ ----------------------------------------------------------------
void JFDMpi2DReflectHighAngle::Calc_Dy(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Dy(j,i)=O1C[0]*WrapAngle(FVal(j-4,i),_F(j,i))+
        O1C[1]*WrapAngle(FVal(j-3,i),_F(j,i))+
        O1C[2]*WrapAngle(FVal(j-2,i),_F(j,i))+
        O1C[3]*WrapAngle(FVal(j-1,i),_F(j,i))+
        O1C[5]*WrapAngle(FVal(j+1,i),_F(j,i))+
        O1C[6]*WrapAngle(FVal(j+2,i),_F(j,i))+
        O1C[7]*WrapAngle(FVal(j+3,i),_F(j,i))+
        O1C[8]*WrapAngle(FVal(j+4,i),_F(j,i));
      _Dy(j,i)/=_dy;
    }
}

// @@ ----------------------------------------------------------------
void JFDMpi2DReflectHighAngle::Calc_Dxx(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Dxx(j,i)=0.0;
      for (int n=0; n<9; n++)
        _Dxx(j,i)+=O2C[n]*WrapAngle(FVal(j,i-4+n),_F(j,i));

      _Dxx(j,i)/=_dxdx;
    }
}

// @@ ----------------------------------------------------------------
void JFDMpi2DReflectHighAngle::Calc_Dyy(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Dyy(j,i)=0.0;
      for (int n=0; n<9; n++)
        _Dyy(j,i)+=O2C[n]*WrapAngle(FVal(j-4+n,i),_F(j,i));

      _Dyy(j,i)/=_dydy;
    }
}
