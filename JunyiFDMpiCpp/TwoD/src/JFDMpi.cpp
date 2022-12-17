#include "JFDMpi.h"

// ## -----------------------------------------------------------
// @@ -----------------------------------------------
JFDMpi2DBase::JFDMpi2DBase(JMpi inJMpi, JMat &in_F) :
  _MpiObj(inJMpi), _F(in_F), _NY(_F.SizeY()), _NX(_F.SizeX()), _Ny(inJMpi.NYLo()),
  _Dx(_NY, _NX), _Dxx(_NY, _NX), _Dy(_NY, _NX), _Dyy(_NY, _NX), _D2(_NY, _NX),
  _FTop(1, _NX), _FBot(1, _NX), _FTop2(1, _NX), _FBot2(1, _NX),
  _dx(_MpiObj.dx()), _dy(_MpiObj.dy()),
  _dxdouble(2.0*_MpiObj.dx()), _dydouble(2.0*_MpiObj.dy()),
  _dxdx(_MpiObj.dx()*_MpiObj.dx()), _dydy(_MpiObj.dy()*_MpiObj.dy())
{}

// @@ -----------------------------------------------
JFDMpi2DBase & JFDMpi2DBase::operator= (const JFDMpi2DBase &in1){

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

  return *this;
}

// @@ -----------------------------------------------
void JFDMpi2DBase::SetF(JMat &MatIn)
{
  _F=MatIn;
}

// @@ -----------------------------------------------
void JFDMpi2DBase::Transfer()
{
  if (_MpiObj.Nnode()==0){
    NDu1=_MpiObj.NYLo()-1;
    NDu2=0;

    for (int i=0; i<_NX; i++){
      _FTop(0,i)=_F(NDu1,i);
      _FBot(0,i)=_F(NDu2,i);
    }

    MPI_Send(_FTop.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()+1,2,MPI_COMM_WORLD);
    MPI_Recv(_FTop.Pointer(),_NX,MPI_DOUBLE,_MpiObj.NLast(),2,MPI_COMM_WORLD,&_status);

    MPI_Send(_FBot.Pointer(),_NX,MPI_DOUBLE,_MpiObj.NLast(),3,MPI_COMM_WORLD);
    MPI_Recv(_FBot.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()+1,3,MPI_COMM_WORLD,&_status);

  } else if (_MpiObj.Nnode()==_MpiObj.NLast()){

    NDu1=_MpiObj.NYLo()-1;
    NDu2=0;

    for (int i=0; i<_NX; i++){
      _FTop(0,i)=_F(NDu1,i);
      _FBot(0,i)=_F(NDu2,i);
    }

    MPI_Recv(_FTop.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()-1,2,MPI_COMM_WORLD,&_status);
    MPI_Send(_FTop.Pointer(),_NX,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
    MPI_Recv(_FBot.Pointer(),_NX,MPI_DOUBLE,0,3,MPI_COMM_WORLD,&_status);
    MPI_Send(_FBot.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()-1,3,MPI_COMM_WORLD);

  } else{
    NDu1=_MpiObj.NYLo()-1;
    NDu2=0;

    for (int i=0; i<_NX; i++){
      _FTop(0,i)=_F(NDu1,i);
      _FBot(0,i)=_F(NDu2,i);
    }

    MPI_Recv(_FTop.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()-1,2,MPI_COMM_WORLD,&_status);
    MPI_Send(_FTop.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()+1,2,MPI_COMM_WORLD);

    MPI_Recv(_FBot.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()+1,3,MPI_COMM_WORLD,&_status);
    MPI_Send(_FBot.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()-1,3,MPI_COMM_WORLD);




  }
}

// @@ -----------------------------------------------
void JFDMpi2DBase::Transfer2()
{
  if (_MpiObj.Nnode()==0){
    NDu1=_MpiObj.NYLo()-2;
    NDu2=1;

    for (int i=0; i<_NX; i++){
      _FTop2(0,i)=_F(NDu1,i);
      _FBot2(0,i)=_F(NDu2,i);
    }

    MPI_Send(_FTop2.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()+1,2,MPI_COMM_WORLD);
    MPI_Recv(_FTop2.Pointer(),_NX,MPI_DOUBLE,_MpiObj.NLast(),2,MPI_COMM_WORLD,&_status);

    MPI_Send(_FBot2.Pointer(),_NX,MPI_DOUBLE,_MpiObj.NLast(),3,MPI_COMM_WORLD);
    MPI_Recv(_FBot2.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()+1,3,MPI_COMM_WORLD,&_status);

  } else if (_MpiObj.Nnode()==_MpiObj.NLast()){

    NDu1=_MpiObj.NYLo()-2;
    NDu2=1;

    for (int i=0; i<_NX; i++){
      _FTop2(0,i)=_F(NDu1,i);
      _FBot2(0,i)=_F(NDu2,i);
    }

    MPI_Recv(_FTop2.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()-1,2,MPI_COMM_WORLD,&_status);
    MPI_Send(_FTop2.Pointer(),_NX,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
    MPI_Recv(_FBot2.Pointer(),_NX,MPI_DOUBLE,0,3,MPI_COMM_WORLD,&_status);
    MPI_Send(_FBot2.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()-1,3,MPI_COMM_WORLD);

  } else{
    NDu1=_MpiObj.NYLo()-2;
    NDu2=1;

    for (int i=0; i<_NX; i++){
      _FTop2(0,i)=_F(NDu1,i);
      _FBot2(0,i)=_F(NDu2,i);
    }

    MPI_Recv(_FTop2.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()-1,2,MPI_COMM_WORLD,&_status);
    MPI_Send(_FTop2.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()+1,2,MPI_COMM_WORLD);
    MPI_Recv(_FBot2.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()+1,3,MPI_COMM_WORLD,&_status);
    MPI_Send(_FBot2.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()-1,3,MPI_COMM_WORLD);

  }
}

// @@ -----------------------------------------------
void JFDMpi2DBase::Calc_Dx()
{
  for (int j=1; j<_MpiObj.NYLo()-1; j++){
		for (int i=1; i<_MpiObj.NX()-1; i++)
			_Dx(j,i)=(_F(j,i+1)-_F(j,i-1))/_dxdouble;
	}
}

// @@ -----------------------------------------------
void JFDMpi2DBase::Calc_Dy()
{
  for (int j=1; j<(_MpiObj.NYLo()-1); j++)
		for (int i=1; i<_NX-1; i++)
			_Dy(j,i)=(_F(j+1,i)-_F(j-1,i))/_dydouble;

}

// @@ -----------------------------------------------
void JFDMpi2DBase::Calc_Dxx()
{
	for (int j=1; j<_MpiObj.NYLo()-1; j++){
		for (int i=1; i<_NX-1; i++)
			_Dxx(j,i)=(_F(j,i+1)+_F(j,i-1)-2.0*_F(j,i))/_dxdx;
	}
}

// @@ -----------------------------------------------
void JFDMpi2DBase::Calc_Dyy()
{
	for (int j=1; j<(_MpiObj.NYLo()-1); j++)
		for (int i=1; i<_NX-1; i++)
			_Dyy(j,i)=(_F(j+1,i)+_F(j-1,i)-2.0*_F(j,i))/_dydy;

}

// @@ -----------------------------------------------
void JFDMpi2DBase::Calc_D2()
{
	for (int j=0; j<_MpiObj.NYLo(); j++)
		for (int i=0; i<_NX; i++)
			_D2(j,i)=_Dxx(j,i)+_Dyy(j,i);
}

// @@ -----------------------------------------------
void JFDMpi2DBase::Calc_O1()
{
	Calc_Dx();
	Calc_Dy();
}

// @@ -----------------------------------------------
void JFDMpi2DBase::Calc_O2ExceptD2()
{
	Calc_Dxx();
	Calc_Dyy();
}

// @@ -----------------------------------------------
void JFDMpi2DBase::Calc_O2()
{
	Calc_Dxx();
	Calc_Dyy();
  Calc_D2();
}

// @@ -----------------------------------------------
void JFDMpi2DBase::Calc_AllExceptD2()
{
  Calc_O1();
  Calc_O2ExceptD2();
}

// @@ -----------------------------------------------
void JFDMpi2DBase::Calc_All()
{
  Calc_O1();
  Calc_O2();
}

// ## =====================================================================
// @@ -----------------------------------------------
JFDMpi2DPeriodicHigh::JFDMpi2DPeriodicHigh(JMpi inJMpi, JMat &in_F) :
	JFDMpi2DBase(inJMpi, in_F),
  _FTop3(1, _NX), _FBot3(1, _NX),
  _FTop4(1, _NX), _FBot4(1, _NX), _Dxy(_NY,_NX) {}

// @@ -----------------------------------------------
void JFDMpi2DPeriodicHigh::Transfer3()
{
  if (_MpiObj.Nnode()==0){
    NDu1=_MpiObj.NYLo()-3;
    NDu2=2;

    for (int i=0; i<_NX; i++){
      _FTop3(0,i)=_F(NDu1,i);
      _FBot3(0,i)=_F(NDu2,i);
    }

    MPI_Send(_FTop3.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()+1,2,MPI_COMM_WORLD);
    MPI_Recv(_FTop3.Pointer(),_NX,MPI_DOUBLE,_MpiObj.NLast(),2,MPI_COMM_WORLD,&_status);
    MPI_Send(_FBot3.Pointer(),_NX,MPI_DOUBLE,_MpiObj.NLast(),3,MPI_COMM_WORLD);
    MPI_Recv(_FBot3.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()+1,3,MPI_COMM_WORLD,&_status);

  } else if (_MpiObj.Nnode()==_MpiObj.NLast()){

    NDu1=_MpiObj.NYLo()-3;
    NDu2=2;

    for (int i=0; i<_NX; i++){
      _FTop3(0,i)=_F(NDu1,i);
      _FBot3(0,i)=_F(NDu2,i);
    }

    MPI_Recv(_FTop3.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()-1,2,MPI_COMM_WORLD,&_status);
    MPI_Send(_FTop3.Pointer(),_NX,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
    MPI_Recv(_FBot3.Pointer(),_NX,MPI_DOUBLE,0,3,MPI_COMM_WORLD,&_status);
    MPI_Send(_FBot3.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()-1,3,MPI_COMM_WORLD);



  } else{
    NDu1=_MpiObj.NYLo()-3;
    NDu2=2;

    for (int i=0; i<_NX; i++){
      _FTop3(0,i)=_F(NDu1,i);
      _FBot3(0,i)=_F(NDu2,i);
    }

    MPI_Recv(_FTop3.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()-1,2,MPI_COMM_WORLD,&_status);
    MPI_Send(_FTop3.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()+1,2,MPI_COMM_WORLD);
    MPI_Recv(_FBot3.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()+1,3,MPI_COMM_WORLD,&_status);
    MPI_Send(_FBot3.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()-1,3,MPI_COMM_WORLD);


  }
}

// @@ -----------------------------------------------
void JFDMpi2DPeriodicHigh::Transfer4()
{
  if (_MpiObj.Nnode()==0){
    NDu1=_MpiObj.NYLo()-4;
    NDu2=3;

    for (int i=0; i<_NX; i++){
      _FTop4(0,i)=_F(NDu1,i);
      _FBot4(0,i)=_F(NDu2,i);
    }

    MPI_Send(_FTop4.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()+1,2,MPI_COMM_WORLD);
    MPI_Recv(_FTop4.Pointer(),_NX,MPI_DOUBLE,_MpiObj.NLast(),2,MPI_COMM_WORLD,&_status);
    MPI_Send(_FBot4.Pointer(),_NX,MPI_DOUBLE,_MpiObj.NLast(),3,MPI_COMM_WORLD);
    MPI_Recv(_FBot4.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()+1,3,MPI_COMM_WORLD,&_status);

  } else if (_MpiObj.Nnode()==_MpiObj.NLast()){

    NDu1=_MpiObj.NYLo()-4;
    NDu2=3;

    for (int i=0; i<_NX; i++){
      _FTop4(0,i)=_F(NDu1,i);
      _FBot4(0,i)=_F(NDu2,i);
    }

    MPI_Recv(_FTop4.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()-1,2,MPI_COMM_WORLD,&_status);
    MPI_Send(_FTop4.Pointer(),_NX,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
    MPI_Recv(_FBot4.Pointer(),_NX,MPI_DOUBLE,0,3,MPI_COMM_WORLD,&_status);
    MPI_Send(_FBot4.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()-1,3,MPI_COMM_WORLD);


  } else{
    NDu1=_MpiObj.NYLo()-4;
    NDu2=3;

    for (int i=0; i<_NX; i++){
      _FTop4(0,i)=_F(NDu1,i);
      _FBot4(0,i)=_F(NDu2,i);
    }

    MPI_Recv(_FTop4.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()-1,2,MPI_COMM_WORLD,&_status);
    MPI_Send(_FTop4.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()+1,2,MPI_COMM_WORLD);
    MPI_Recv(_FBot4.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()+1,3,MPI_COMM_WORLD,&_status);
    MPI_Send(_FBot4.Pointer(),_NX,MPI_DOUBLE,_MpiObj.Nnode()-1,3,MPI_COMM_WORLD);


  }
}

// @@ -----------------------------------------------
void JFDMpi2DPeriodicHigh::TransferAll(){

  Transfer();
  Transfer2();
  Transfer3();
  Transfer4();
}

// @@ ------------------------------------------------------------------
double JFDMpi2DPeriodicHigh::FVal(const int &y, const int &x){


  if ((x>=0)&&(x<_NX))
    xTemp=x;
  else if (x<0)
    xTemp=x+_NX;
  else
    xTemp=x-_NX;

  if ((y>=0)&&(y<_Ny))
    return _F(y,xTemp);
  else if (y==-1)
    return _FTop(0,xTemp);
  else if (y==-2)
    return _FTop2(0,xTemp);
  else if (y==-3)
    return _FTop3(0,xTemp);
  else if (y==-4)
    return _FTop4(0,xTemp);
  else if (y==_Ny)
    return _FBot(0,xTemp);
  else if (y==_Ny+1)
    return _FBot2(0,xTemp);
  else if (y==_Ny+2)
    return _FBot3(0,xTemp);
  else if (y==_Ny+3)
    return _FBot4(0,xTemp);
  else
    return _F(y,xTemp);

	// return _F(yTemp,xTemp);
}

// @@ ------------------------------------------------------------------
double JFDMpi2DPeriodicHigh::DyVal(const int &y, const int &x){

  if ((x>=0)&&(x<_NX))
    xTemp=x;
  else if (x<0)
    xTemp=x+_NX;
  else
    xTemp=x-_NX;

	return _Dy(y,xTemp);
}

// @@ ----------------------------------------------------------------
void JFDMpi2DPeriodicHigh::Calc_Dx(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Dx(j,i)=O1C[0]*FVal(j,i-4)+
        O1C[1]*FVal(j,i-3)+
        O1C[2]*FVal(j,i-2)+
        O1C[3]*FVal(j,i-1)+
        O1C[5]*FVal(j,i+1)+
        O1C[6]*FVal(j,i+2)+
        O1C[7]*FVal(j,i+3)+
        O1C[8]*FVal(j,i+4);
      _Dx(j,i)/=_dx;
    }
}

// @@ ----------------------------------------------------------------
void JFDMpi2DPeriodicHigh::Calc_Dy(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Dy(j,i)=O1C[0]*FVal(j-4,i)+
        O1C[1]*FVal(j-3,i)+
        O1C[2]*FVal(j-2,i)+
        O1C[3]*FVal(j-1,i)+
        O1C[5]*FVal(j+1,i)+
        O1C[6]*FVal(j+2,i)+
        O1C[7]*FVal(j+3,i)+
        O1C[8]*FVal(j+4,i);
      _Dy(j,i)/=_dy;
    }
}

// @@ ----------------------------------------------------------------
void JFDMpi2DPeriodicHigh::Calc_Dxy(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Dxy(j,i)=O1C[0]*DyVal(j,i-4)+
        O1C[1]*DyVal(j,i-3)+
        O1C[2]*DyVal(j,i-2)+
        O1C[3]*DyVal(j,i-1)+
        O1C[5]*DyVal(j,i+1)+
        O1C[6]*DyVal(j,i+2)+
        O1C[7]*DyVal(j,i+3)+
        O1C[8]*DyVal(j,i+4);
      _Dxy(j,i)/=_dx;
    }
}

// @@ ----------------------------------------------------------------
void JFDMpi2DPeriodicHigh::Calc_Dxx(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Dxx(j,i)=0.0;
      for (int n=0; n<9; n++)
        _Dxx(j,i)+=O2C[n]*FVal(j,i-4+n);

      _Dxx(j,i)/=_dxdx;
    }
}

// @@ ----------------------------------------------------------------
void JFDMpi2DPeriodicHigh::Calc_Dyy(){
  for (int j=0; j<_Ny; j++)
    for (int i=0; i<_NX; i++){
      _Dyy(j,i)=0.0;
      for (int n=0; n<9; n++)
        _Dyy(j,i)+=O2C[n]*FVal(j-4+n,i);

      _Dyy(j,i)/=_dydy;
    }
}

// @@ ----------------------------------------------------------------
void JFDMpi2DPeriodicHigh::Calc_O1(){
  Calc_Dx();
  Calc_Dy();
  Calc_Dxy();
}

// @@ ----------------------------------------------------------------
void JFDMpi2DPeriodicHigh::Calc_O1ExceptXY(){
  Calc_Dx();
  Calc_Dy();
}

// @@ ----------------------------------------------------------------
void JFDMpi2DPeriodicHigh::Calc_O2ExceptD2(){
  Calc_Dxx();
  Calc_Dyy();
}

// @@ ----------------------------------------------------------------
void JFDMpi2DPeriodicHigh::Calc_O2(){
  Calc_Dxx();
  Calc_Dyy();
  Calc_D2();
}

// @@ ----------------------------------------------------------------
void JFDMpi2DPeriodicHigh::Calc_AllExceptD2(){
  Calc_O1();
  Calc_O2ExceptD2();
}

// @@ ----------------------------------------------------------------
void JFDMpi2DPeriodicHigh::Calc_AllExceptXY(){
  Calc_O1ExceptXY();
  Calc_O2();
}

// @@ ----------------------------------------------------------------
void JFDMpi2DPeriodicHigh::Calc_All(){
  Calc_O1();
  Calc_O2();
}

// @@ -----------------------------------------------
JFDMpi2DPeriodicHigh & JFDMpi2DPeriodicHigh::operator= (const JFDMpi2DPeriodicHigh &in1){

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

// ## =====================================================================
// @@ -----------------------------------------------
JFDMpi2DReflectHigh::JFDMpi2DReflectHigh(JMpi inJMpi, JMat &in_F) :
	JFDMpi2DPeriodicHigh(inJMpi, in_F) {}

// @@ -----------------------------------------------
JFDMpi2DReflectHigh & JFDMpi2DReflectHigh::operator= (const JFDMpi2DReflectHigh &in1){

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

// @@ ------------------------------------------------------------------
double JFDMpi2DReflectHigh::FVal(const int &y, const int &x){

//===================================================
  if ((_MpiObj.Nnode()!=0)&&(_MpiObj.Nnode()!=_MpiObj.NLast())){

    if ((x>=0)&&(x<_NX))
      xTemp=x;
    else if (x<0)
      xTemp=-x;
    else
      xTemp=_NX+_NX-2-x;

    if ((y>=0)&&(y<_Ny))
      return _F(y,xTemp);
    else if (y==-1)
      return _FTop(0,xTemp);
    else if (y==-2)
      return _FTop2(0,xTemp);
    else if (y==-3)
      return _FTop3(0,xTemp);
    else if (y==-4)
      return _FTop4(0,xTemp);
    else if (y==_Ny)
      return _FBot(0,xTemp);
    else if (y==_Ny+1)
      return _FBot2(0,xTemp);
    else if (y==_Ny+2)
      return _FBot3(0,xTemp);
    else if (y==_Ny+3)
      return _FBot4(0,xTemp);
    else
      return _F(y,xTemp);

  } else if (_MpiObj.Nnode()==0) {

    if ((x>=0)&&(x<_NX))
      xTemp=x;
    else if (x<0)
      xTemp=-x;
    else
      xTemp=_NX+_NX-2-x;

    if ((y>=0)&&(y<_Ny))
      return _F(y,xTemp);
    else if (y==-1)
      return _F(1,xTemp);
    else if (y==-2)
      return _F(2,xTemp);
    else if (y==-3)
      return _F(3,xTemp);
    else if (y==-4)
      return _F(4,xTemp);
    else if (y==_Ny)
      return _FBot(0,xTemp);
    else if (y==_Ny+1)
      return _FBot2(0,xTemp);
    else if (y==_Ny+2)
      return _FBot3(0,xTemp);
    else if (y==_Ny+3)
      return _FBot4(0,xTemp);
    else
      return _F(y,xTemp);

  } else {

    if ((x>=0)&&(x<_NX))
      xTemp=x;
    else if (x<0)
      xTemp=-x;
    else
      xTemp=_NX+_NX-2-x;

    if ((y>=0)&&(y<_Ny))
      return _F(y,xTemp);
    else if (y==-1)
      return _FTop(0,xTemp);
    else if (y==-2)
      return _FTop2(0,xTemp);
    else if (y==-3)
      return _FTop3(0,xTemp);
    else if (y==-4)
      return _FTop4(0,xTemp);
    else if (y==_Ny)
      return _F(_Ny-2,xTemp);
    else if (y==_Ny+1)
      return _F(_Ny-3,xTemp);
    else if (y==_Ny+2)
      return _F(_Ny-4,xTemp);
    else if (y==_Ny+3)
      return _F(_Ny-5,xTemp);
    else
      return _F(y,xTemp);

  }

	// return _F(yTemp,xTemp);
}

// @@ ------------------------------------------------------------------
double JFDMpi2DReflectHigh::DyVal(const int &y, const int &x){

  if ((x>=0)&&(x<_NX))
    xTemp=x;
  else if (x<0)
    xTemp=-x;
  else
    xTemp=_NX+_NX-2-x;

	return _Dy(y,xTemp);
}
