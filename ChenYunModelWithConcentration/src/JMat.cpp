#include "JMat.h"

// ## JMat-----------------------------------------------------------------
// @@ ---------------------------------------------------------------
JMat::JMat()
{
  _NY=1;
  _NX=1;
  _NXY=_NY*_NX;
  _XPtr=new double[_NXY];

  for (int j=0; j<_NY; j++)
    for (int i=0; i<_NX; i++)
      *(_XPtr+i+j*_NX)=0.0;
}

// @@ ---------------------------------------------------------------
JMat::JMat(const JMat &inJMat) :
  _NY(inJMat._NY), _NX(inJMat._NX), _NXY(inJMat._NXY)
{
	_XPtr=new double[_NXY];

	for (int j=0; j<_NY; j++)
		for (int i=0; i<_NX; i++)
			*(_XPtr+i+j*_NX)=*(inJMat._XPtr+i+j*_NX);
}

// @@ ---------------------------------------------------------------
JMat::JMat(const int &in_NY, const int &in_NX) :
  _NY(in_NY), _NX(in_NX), _NXY(_NY*_NX)
{
  _XPtr=new double[_NXY];

  for (int j=0; j<_NY; j++)
    for (int i=0; i<_NX; i++)
      *(_XPtr+i+j*_NX)=0.0;
}

// @@ ---------------------------------------------------------------
JMat::~JMat()
{
	delete[] _XPtr;
}

// @@ ---------------------------------------------------------------
JMat & JMat::operator= (const JMat & inJMat)
{
  if (inJMat._NXY!=_NXY){
		_NX=inJMat._NX;
		_NY=inJMat._NY;
		_NXY=inJMat._NXY;

		delete [] _XPtr;
		_XPtr=new double[_NXY];
		std::cout << "Warning Pointer Change" <<std::endl;
	}

	for (int j=0; j<_NY; j++)
		for (int i=0; i<_NX; i++)
			*(_XPtr+i+j*_NX)=*(inJMat._XPtr+i+j*_NX);

	return *this;
}

// @@ ---------------------------------------------------------------
double & JMat::operator() (const int &iny, const int &inx)
{
  return *(_XPtr+inx+iny*_NX);
}

// @@ ---------------------------------------------------------------
void JMat::Print()
{
  for (int j=0; j<_NY; j++){
		for (int i=0; i<_NX-1; i++)
			std::cout<<*(_XPtr+i+j*_NX)<<',';
		std::cout<<*(_XPtr+(_NX-1)+j*_NX)<<std::endl;
	}
}

// @@ ---------------------------------------------------------------
void JMat::Resize(const int &iny,const  int &inx){
	_NX=inx;
	_NY=iny;
	_NXY=_NY*_NX;
	delete [] _XPtr;
	_XPtr=new double[_NXY];

	for (int j=0; j<_NY; j++)
		for (int i=0; i<_NX; i++)
			*(_XPtr+i+j*_NX)=0.0;
}

// @@ ---------------------------------------------------------------
void JMat::Replace(const int &iny, const int &inx, const double &inval)
{
  *(_XPtr+inx+iny*_NX)=inval;
}

// @@ ---------------------------------------------------------------
void JMat::Zero()
{
  for (int j=0;j<_NY;j++)
    for (int i=0;i<_NX;i++)
      Replace(j, i, 0.0);
}
