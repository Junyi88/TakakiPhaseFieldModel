#ifndef JMATCLASS_H
#define JMATCLASS_H

#include <iostream>
#include <string>
#include <istream>
#include <fstream>
#include <sstream>
#include <cmath>

// ## -----------------------------------------------------------------
class JMat{
public:
  // Constructors, Destructors and Key operators
  JMat();
  JMat(const JMat &inJMat);
  JMat(const int &in_NY, const int &in_NX);
  ~JMat();

  JMat & operator= (const JMat & inJMat); //Write to operator
	double & operator() (const int &iny, const int &inx); // Output operator

  // Member functions
  void Print();
  void Resize(const int &iny,const  int &inx);
  void Replace(const int &iny, const int &inx, const double &inval);
  void Zero();

  int SizeX(){return _NX;};
  int SizeY(){return _NY;};
  int Size(){return _NXY;};

  double Value(const int &iny, const int &inx){return *(_XPtr+inx+iny*_NX);};

  double * Pointer(){return _XPtr;};
  
// Members
protected:
  int _NY,_NX,_NXY;
	double * _XPtr;
};


#endif
