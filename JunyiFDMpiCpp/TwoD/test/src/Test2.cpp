#include<iostream>
#include<mpi.h>
#include<cmath>
#include <string>
#include <vector>
#include <istream>
#include <fstream>
#include <sstream>
#include "JMat.h"
#include "JMpi.h"
#include "JFDMpi.h"
#include "AuxFunctions.h"

int main(int argc, char ** argv){

  //================================================================
  // # Initialise the MPI
  int NPrs, Nnode;
  MPI_Status status;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &NPrs);
  MPI_Comm_rank(MPI_COMM_WORLD, &Nnode);

  std::string BuffString1;
  std::ifstream BuffInputStream1;
  int BuffInt;

  std::string HeaderName;
  std::ofstream LogFile, BuffOutputFile1;
  std::string filename;

  std::vector<std::string> InputFilesVector, InitialFilesVector;
  int NY=201;
  int NX=201;
  double dy=1.0;
  double dx=1.0;

  JMpi MpiObj(Nnode, NPrs, NY, NX, dy, dx);
  JMat FullBuffMat(MpiObj.NYFu(), MpiObj.NX()), F(MpiObj.NYGl(), MpiObj.NX());



  //===============================================================
  HeaderName=argv[1];
  filename=HeaderName + "_Node" + std::to_string(Nnode) + ".log";
  LogFile.open(filename);
  LogFile << "Open File Node " << Nnode << std::endl;
  LogFile << "NYLo = " << MpiObj.NYLo() << std::endl;
  LogFile << "NYGl = " << MpiObj.NYGl() << std::endl;
  LogFile << "NYShort = " << MpiObj.NYShort()<< std::endl;
  LogFile << "YStart = " << MpiObj.YStart()<< std::endl;
  LogFile << "PosShortLast = " << MpiObj.PosShortLast()<< std::endl;
  //=====================================================================
  LogFile << "================================="  << std::endl;

  filename=argv[2];
  LogFile << "Reading Input File = " << filename << std::endl;

  ReadTextFile(FullBuffMat.Pointer(), filename, MpiObj.NX(), MpiObj.NYFu(), ',');
  Splitter(F, FullBuffMat, MpiObj);
  LogFile << "=====Build MyFD"  << std::endl;
  JFDMpi2DReflectHigh MyFD(MpiObj, F);
  LogFile << "=====Build MyFD Done"  << std::endl;
  MyFD.TransferAll();
  LogFile << "=====Transfer"  << std::endl;
  MyFD.Calc_All();
  LogFile << "=====Calculation"  << std::endl;

  filename=HeaderName + "_F.csv";
  WriteMPITextFile(MyFD.FP(), filename, MpiObj);

  filename=HeaderName + "_Dx.csv";
  WriteMPITextFile(MyFD.DxP(), filename, MpiObj);
  filename=HeaderName + "_Dy.csv";
  WriteMPITextFile(MyFD.DyP(), filename, MpiObj);

  filename=HeaderName + "_Dxx.csv";
  WriteMPITextFile(MyFD.DxxP(), filename, MpiObj);
  filename=HeaderName + "_Dyy.csv";
  WriteMPITextFile(MyFD.DyyP(), filename, MpiObj);

  filename=HeaderName + "_Dxy.csv";
  WriteMPITextFile(MyFD.DxyP(), filename, MpiObj);
  filename=HeaderName + "_D2.csv";
  WriteMPITextFile(MyFD.D2P(), filename, MpiObj);

  //===============================================================
  LogFile << "Finished Run \n" << std::endl;
  LogFile.close();
  MPI_Finalize();
  return 0;
}
