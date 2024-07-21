#ifndef JFDMPIANGLE_H
#define JFDMPIANGLE_H

#define _USE_MATH_DEFINES // For VISUAL STUDIO
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


#include <iostream>
#include <string>
#include <istream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <mpi.h>
#include "JMpi.h"
#include "JMat.h"
#include "JFDMpi.h"

// ## -----------------------------------------------------------------
class JFDMpi2DReflectHighAngle : public JFDMpi2DReflectHigh {

public:
	//JFDMpi2DPeriodicHigh(){ std::cout<<"JFDMpi2DPeriodicHigh empty construction error"<<std::endl; };
	JFDMpi2DReflectHighAngle(JMpi inJMpi, JMat &in_F);
	JFDMpi2DReflectHighAngle & operator= (const JFDMpi2DReflectHighAngle &in1); //Write to operator

	//--------------------------------------------------------------
	void Calc_Dx() override;
  void Calc_Dy() override;
	void Calc_Dxx() override;
  void Calc_Dyy() override;

	// double FVal(const int &y, const int &x) override;
	double WrapAngle(const double &Next, const double &Ref);

	void CleanUpAllAngles();

private:
	const double Pi2=2*M_PI;
	double ReferenceAngle;
	double AngleBuffer;
};

#endif
