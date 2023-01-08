#ifndef JFDMPI_H
#define JFDMPI_H

#include <iostream>
#include <string>
#include <istream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <mpi.h>
#include "JMpi.h"
#include "JMat.h"

// ## -----------------------------------------------------------------
class JFDMpi2DBase{
public:
	//JFDMpi2DBase(){ std::cout<<"JFDMpi2DBase empty construction error"<<std::endl; };
	JFDMpi2DBase(JMpi inJMpi, JMat &in_F);
	JFDMpi2DBase & operator= (const JFDMpi2DBase &in1); //Write to operator

// -----------------------------------
// Key Functions
public:
	void SetF(JMat &MatIn);

  void Transfer();
	void Transfer2();
	//
  virtual void Calc_Dx();
  virtual void Calc_Dy();

  virtual void Calc_Dxx();
  virtual void Calc_Dyy();
	void Calc_D2();

  virtual void Calc_O1();
  virtual void Calc_O2ExceptD2();
	virtual void Calc_O2();

  virtual void Calc_AllExceptD2();
  virtual void Calc_All();

// Getter Functions
	double Dx(const int &y, const int &x){return _Dx.Value(y,x);};
	double Dy(const int &y, const int &x){return _Dy.Value(y,x);};

	double Dxx(const int &y, const int &x){return _Dxx.Value(y,x);};
	double Dyy(const int &y, const int &x){return _Dyy.Value(y,x);};
	double D2(const int &y, const int &x){return _D2.Value(y,x);};

	JMat * FP(){return &_F;};
	JMat * DxP(){return &_Dx;};
	JMat * DyP(){return &_Dy;};
	JMat * DxxP(){return &_Dxx;};
	JMat * DyyP(){return &_Dyy;};
	JMat * D2P(){return &_D2;};

	JMat & FR(){return _F;};
	JMat & DxR(){return _Dx;};
	JMat & DyR(){return _Dy;};
	JMat & DxxR(){return _Dxx;};
	JMat & DyyR(){return _Dyy;};
	JMat & D2R(){return _D2;};

	double * DxPointer(){return _Dx.Pointer();};
	double * DyPointer(){return _Dy.Pointer();};
	double * DxxPointer(){return _Dxx.Pointer();};
	double * DyyPointer(){return _Dyy.Pointer();};
	double * D2Pointer(){return _D2.Pointer();};
// -----------------------------------
// Members
protected:
  JMpi _MpiObj;
  JMat &_F;

	int _NY, _NX, _Ny;
  JMat _Dx, _Dxx, _Dy, _Dyy, _D2;
  JMat _FTop, _FBot, _FTop2, _FBot2;
  JMat _FTopSend, _FBotSend, _FTop2Send, _FBot2Send;
  double _dx,_dy,_dxdouble, _dydouble, _dxdx, _dydy;

  MPI_Status _status;

	int NDu1, NDu2, NDu3, NDu4, NDu5;
private:

};

// ## -----------------------------------------------------------------
class JFDMpi2DPeriodicHigh : public JFDMpi2DBase {

public:
	//JFDMpi2DPeriodicHigh(){ std::cout<<"JFDMpi2DPeriodicHigh empty construction error"<<std::endl; };
	JFDMpi2DPeriodicHigh(JMpi inJMpi, JMat &in_F);
	JFDMpi2DPeriodicHigh & operator= (const JFDMpi2DPeriodicHigh &in1); //Write to operator


	// Key Functions
	void Calc_Dx() override;
  void Calc_Dy() override;
	void Calc_Dxy();

  void Calc_Dxx() override;
  void Calc_Dyy() override;

  void Calc_O1() override;
	void Calc_O1ExceptXY();
  void Calc_O2ExceptD2() override;
	void Calc_O2() override;

  void Calc_AllExceptD2() override;
	void Calc_AllExceptXY();
  void Calc_All() override;

	void Transfer3();
	void Transfer4();
	void TransferAll();

	virtual double FVal(const int &y, const int &x);
	virtual double DyVal(const int &y, const int &x);

	// Getter Functions
	double Dxy(const int &y, const int &x){return _Dxy.Value(y,x);};
	JMat * DxyP(){return &_Dxy;};
	double * DxyPointer(){return _Dxy.Pointer();};
protected:
	const double O1C[9] = { 1.0/280.0, -4.0/105.0, 1.0/5.0, -4.0/5.0, 0.0,
	  	4.0/5.0, -1.0/5.0, 4.0/105.0, -1.0/280.0};
	const double O2C[9] = {-1.0/560.0, 8.0/315.0, -1.0/5.0, 8.0/5.0, -205.0/72.0,
	   	8.0/5.0, -1.0/5.0, 8.0/315.0, -1.0/560.0};
	JMat _FTop3, _FBot3, _FTop4, _FBot4;
	JMat _FTop3Send, _FBot3Send, _FTop4Send, _FBot4Send;

	JMat _Dxy;
	int xTemp, yTemp;
private:

};

// ## -----------------------------------------------------------------
class JFDMpi2DReflectHigh : public JFDMpi2DPeriodicHigh {

public:
	//JFDMpi2DPeriodicHigh(){ std::cout<<"JFDMpi2DPeriodicHigh empty construction error"<<std::endl; };
	JFDMpi2DReflectHigh(JMpi inJMpi, JMat &in_F);
	JFDMpi2DReflectHigh & operator= (const JFDMpi2DReflectHigh &in1); //Write to operator

	double FVal(const int &y, const int &x) override;
	double DyVal(const int &y, const int &x) override;

};


#endif
