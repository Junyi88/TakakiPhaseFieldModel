#ifndef JMPI_H
#define JMPI_H

#include <iostream>
#include <string>
#include <istream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <mpi.h>

// ## -----------------------------------------------------------------
class JMpi{
public:
  // Constructors, Destructors and Key operators
  JMpi();
  JMpi(const JMpi &inJMpi);
  JMpi(const int & inNnode, const int & inNPrs, const int & inNYFu, const int &inNX);
  JMpi(const int & inNnode, const int & inNPrs, const int & inNYFu, const int &inNX,
    const double &indy, const double &indx);
  JMpi(const int & inNYFu, const int &inNX, const double &indy, const double &indx);

  JMpi & operator= (const JMpi &inJMpi);

  // Useful functions
  void Set(const int & inNnode, const int & inNPrs, const int & inNYFu, const int &inNX);
  void wait();
  void CalcSizes();

  // Access to protected parts
  int Nnode(){return _Nnode;};
  int NPrs(){return _NPrs;};
  int Nwork(){return _Nwork;};
  int NLast(){return _NLast;};
  double dy(){return _dy;};
  double dx(){return _dx;};

  int NX(){return _NX;};
  int NYFu(){return _NYFu;};
  int NYLo(){return _NYLo;};
  int NYGl(){return _NYGl;};
  int NYShort(){return _NYShort;};
  int NXYGl(){return _NXYGl;};
  int PosShortLast(){return _PosShortLast;};
  int YStart(){return _YStart;};

  MPI_Status * MpiStatusPointer(){return &(_status);};
  int * WaitFlagPointer(){return &(_WaitFlag);};

protected:
  int _Nnode, _NPrs, _Nwork;
  int _NYFu, _NX, _NXYFu;
  int _NYGl, _NXYGl;
  int _NYLo, _NXYLo;
  int _NYLa, _NXYLa;
  int _WaitFlag;
  int _NLast;
  double _dy, _dx;
  MPI_Status _status;


  int _NYShort, _NXYShort;
  int _PosShortLast;
  int _YStart;
};

#endif
