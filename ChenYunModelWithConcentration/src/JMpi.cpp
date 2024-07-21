#include "JMpi.h"

// ## JMpi-----------------------------------------------------------------
// @@ ---------------------------------------------------------------
JMpi::JMpi() :
  _Nnode(0), _NPrs(0), _Nwork(0), _NYFu(0), _NX(0), _NXYFu(0),
  _NYGl(0) , _NXYGl(0), _NYLo(0), _NXYLo(0), _NYLa(0), _NXYLa(0),
  _WaitFlag(0), _NLast(0), _dy(1.0), _dx(1.0)
{
  MPI_Comm_size(MPI_COMM_WORLD, &_NPrs);
  MPI_Comm_rank(MPI_COMM_WORLD, &_Nnode);
}

// @@ ---------------------------------------------------------------
JMpi::JMpi(const JMpi &inJMpi) :
  _Nnode(inJMpi._Nnode), _NPrs(inJMpi._NPrs), _Nwork(inJMpi._Nwork),
  _NYFu(inJMpi._NYFu), _NX(inJMpi._NX), _NXYFu(inJMpi._NXYFu),
  _NYGl(inJMpi._NYGl) , _NXYGl(inJMpi._NXYGl), _NYLo(inJMpi._NYLo),
  _NXYLo(inJMpi._NXYLo), _NYLa(inJMpi._NYLa), _NXYLa(inJMpi._NXYLa),
  _WaitFlag(inJMpi._WaitFlag), _NLast(inJMpi._NLast), _dy(inJMpi._dy),
  _dx(inJMpi._dx)
{}

// @@ ---------------------------------------------------------------
JMpi::JMpi(const int & inNnode, const int & inNPrs, const int & inNYFu, const int &inNX) :
  _Nnode(inNnode), _NPrs(inNPrs), _Nwork(0), _NYFu(inNYFu), _NX(inNX), _NXYFu(0),
  _NYGl(0) , _NXYGl(0), _NYLo(0), _NXYLo(0), _NYLa(0), _NXYLa(0),
  _WaitFlag(0), _NLast(0), _dy(1.0), _dx(1.0)
{
  CalcSizes();
}

// @@ ---------------------------------------------------------------
JMpi::JMpi(const int & inNnode, const int & inNPrs, const int & inNYFu, const int &inNX,
  const double &indy, const double &indx) :
  _Nnode(inNnode), _NPrs(inNPrs), _Nwork(0), _NYFu(inNYFu), _NX(inNX), _NXYFu(0),
  _NYGl(0) , _NXYGl(0), _NYLo(0), _NXYLo(0), _NYLa(0), _NXYLa(0),
  _WaitFlag(0), _NLast(0), _dy(indy), _dx(indx)
{
  CalcSizes();
}

// @@ ---------------------------------------------------------------
JMpi::JMpi(const int & inNYFu, const int &inNX, const double &indy, const double &indx) :
  _Nnode(0), _NPrs(0), _Nwork(0), _NYFu(inNYFu), _NX(inNX), _NXYFu(0),
  _NYGl(0) , _NXYGl(0), _NYLo(0), _NXYLo(0), _NYLa(0), _NXYLa(0),
  _WaitFlag(0), _NLast(0), _dy(indy), _dx(indx)
{
  MPI_Comm_size(MPI_COMM_WORLD, &_NPrs);
  MPI_Comm_rank(MPI_COMM_WORLD, &_Nnode);
  CalcSizes();
}

// @@ ---------------------------------------------------------------
JMpi & JMpi::operator= (const JMpi &inJMpi)
{
  _Nnode=inJMpi._Nnode;
  _NPrs=inJMpi._NPrs;
  _Nwork=inJMpi._Nwork;
  _NYFu=inJMpi._NYFu;
  _NX=inJMpi._NX;
  _NXYFu=inJMpi._NXYFu;
  _NYGl=inJMpi._NYGl;
  _NXYGl=inJMpi._NXYGl;
  _NYLo=inJMpi._NYLo;
  _NXYLo=inJMpi._NXYLo;
  _NYLa=inJMpi._NYLa;
  _NXYLa=inJMpi._NXYLa;
  _WaitFlag=inJMpi._WaitFlag;
  _NLast=inJMpi._NLast;
  _dy=inJMpi._dy;
  _dx=inJMpi._dx;

  _NYShort=inJMpi._NYShort;
  _NXYShort=inJMpi._NXYShort;
  _PosShortLast=inJMpi._PosShortLast;
  _YStart=inJMpi._YStart;
  return *this;
}

// @@ ---------------------------------------------------------------
void JMpi::Set(const int & inNnode, const int & inNPrs, const int & inNYFu, const int &inNX){
  _Nnode=inNnode;
  _NPrs=inNPrs;
  _NYFu=inNYFu;
  _NX=inNX;

  CalcSizes();
}

// @@ ---------------------------------------------------------------
void JMpi::wait(){
  if (_Nnode==0){
    for (int n=1; n<_NPrs; n++){
			MPI_Send(&_WaitFlag, 1, MPI_INT, n, 1001, MPI_COMM_WORLD);
      MPI_Recv(&_WaitFlag, 1, MPI_INT, n, 1001, MPI_COMM_WORLD,&_status);
		}
  }
  else{
		MPI_Send(&_WaitFlag, 1, MPI_INT, 0, 1001, MPI_COMM_WORLD);
		MPI_Recv(&_WaitFlag, 1, MPI_INT, 0, 1001, MPI_COMM_WORLD,&_status);
	}
}


// @@ ---------------------------------------------------------------
void JMpi::CalcSizes()
{
  _Nwork=_NPrs;
  _NLast=_Nwork-1;

  if (_NYFu%_Nwork==0){
    _NYShort=_NYFu/_Nwork;
    _NYGl=_NYShort;
    _NYLa=_NYGl;
  } else{
    _NYShort=_NYFu/_Nwork;
    _NYGl=_NYShort+1;
    _NYLa=_NYFu-(_NLast)*_NYGl;
  }

  if (_NYLa<3)
    if (_Nnode==0)
      std::cout<<"Error NYLa < 3 !" <<std::endl;

  _PosShortLast=_NYFu-(_NPrs*_NYShort); //Dodgy

  if (_Nnode<_PosShortLast){
    _NYLo=_NYGl;
    _YStart=_Nnode*_NYGl;
  } else{
    _NYLo=_NYShort;
    _YStart=_PosShortLast*_NYGl+(_Nnode-_PosShortLast)*_NYShort;
  }




  _NXYFu=_NYFu*_NX;
  _NXYGl=_NYGl*_NX;
  _NXYLo=_NYLo*_NX;
  _NXYLa=_NYLa*_NX;
  _NXYShort=_NYShort*_NX;

  _WaitFlag=0;
}
