//#include "stdafx.h"
#include <iostream>
#include <string>
#include <istream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>

#include "mpi.h"

#include "JMpi.h"
#include "JMat.h"
#include "JFDMpi.h"
#include "JFDMpiAngle.h"
#include "AuxFunctions.h"

#include "TakPhase.h"
#include "TakQuaternion.h"
#include "BasicChemPotential.h"

#include "TakACBulkEnergy.h"
#include "TakACWallEnergy.h"
#include "TakACGradEnergy.h"
#include "TakACTOriEnergyQuaternion.h"
#include "TakACChemEnergy.h"

#include "TakakiSolverQuaternionCon.h"



// #include "GClass.h"

int main(int argc, char ** argv){

  const int NInputParameters=20;
  //================================================================
  // # Initialise the MPI and system
  int NPrs, Nnode;
  MPI_Status status;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &NPrs);
  MPI_Comm_rank(MPI_COMM_WORLD, &Nnode);

  std::string HeaderName, filename;
  std::ofstream LogFile;

  HeaderName=argv[1];
  filename=HeaderName + "_Node" + std::to_string(Nnode) + ".log";
  LogFile.open(filename);
  LogFile << "Open File Node " << Nnode << std::endl;

  // Setup Initial File Name and Inputs
  std::string BufferString;
  std::string InitialConditionFile;
  std::vector <std::string> InitialConditionFileList;
  std::string InputFile;
  std::ifstream BufferInputStream1;

  InputFile=argv[2];


  //====================================================================
  LogFile << "========================================" << std::endl;
  LogFile << "Reading InputFile On Node" << Nnode << std::endl;


  LogFile << "Input File: " << InputFile <<std::endl;
  LogFile << "Reading Initial Files List Completed \n" << std::endl;

  //====================================================================
  LogFile << "========================================" << std::endl;
  int NY, NX, ntStart, ntEnd, WriteCount;
  double dy, dx, dt, MinAngle0, BVec, mu, MTheta0, InvPhiMin, inMPhiConst,
    alpha, Wa, S, T, b, M0, Q, sg;

  NY = 100;
  NX = 100;
  dy = 0.1;
  dx = 0.1;

  //====================================================================
  LogFile << "========================================" << std::endl;
  LogFile << "Setup MPI Object" << std::endl;

  JMpi MPIOBJ(Nnode, NPrs, NY, NX, dy, dx);
  JMat BufferFull(NY, NX), Eta0(MPIOBJ.NYLo(), NX),
      Q10(MPIOBJ.NYLo(), NX),
	  Q20(MPIOBJ.NYLo(), NX),
	  Q30(MPIOBJ.NYLo(), NX),
	  Q40(MPIOBJ.NYLo(), NX),
	  Rho0(MPIOBJ.NYLo(), NX);

  JMat Con0(MPIOBJ.NYLo(), NX);

  LogFile << "Setup MPI Object Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Reading Initial Conditions Eta =  " << InputFile << std::endl;
  BufferString=InputFile;
  ReadTextFile(BufferFull.Pointer() , BufferString, NX, NY, ',');
  Splitter(Eta0, BufferFull,  MPIOBJ);
  LogFile << "Reading Initial Conditions Eta Completed \n" << std::endl;

  // *********************************X****************************************
  // FDClass=JFDMpi2DReflectHigh;
  // FDAngleClass=JFDMpi2DReflectHighAngle;
  LogFile << "========================================" << std::endl;
  LogFile << "Setup Phase Class " << std::endl;
  TakPhase<JFDMpi2DReflectHigh> Phi(MPIOBJ, Eta0);
  LogFile << "Setup Phase Class Completed \n" << std::endl;



  LogFile << "========================================" << std::endl;
  LogFile << "Do 1 Step " << std::endl;
  Phi.Calc_All();
  LogFile << "Done \n" << std::endl;

  //===============================================================
  // Write outputs For First STep
  LogFile << "========================================" << std::endl;
  LogFile << "Write Outputs - Phi" << std::endl;

  BufferString="Phi_0.csv";
  WriteMPITextFile(Phi.FP(), BufferString, MPIOBJ);
  BufferString="Phi2_0.csv";
  WriteMPITextFile(Phi.F2P(), BufferString, MPIOBJ);
  BufferString="Phi3_0.csv";
  WriteMPITextFile(Phi.F3P(), BufferString, MPIOBJ);
  BufferString="Phi4_0.csv";
  WriteMPITextFile(Phi.F4P(), BufferString, MPIOBJ);
  BufferString="Phi5_0.csv";
  WriteMPITextFile(Phi.F5P(), BufferString, MPIOBJ);

  BufferString="PhiDx_0.csv";
  WriteMPITextFile(Phi.DxP(), BufferString, MPIOBJ);
  BufferString="PhiDy_0.csv";
  WriteMPITextFile(Phi.DyP(), BufferString, MPIOBJ);
  BufferString="PhiDxx_0.csv";
  WriteMPITextFile(Phi.DxxP(), BufferString, MPIOBJ);
  BufferString="PhiDyy_0.csv";
  WriteMPITextFile(Phi.DyyP(), BufferString, MPIOBJ);
  BufferString="PhiDxy_0.csv";
  WriteMPITextFile(Phi.DxyP(), BufferString, MPIOBJ);
  BufferString="PhiD2_0.csv";
  WriteMPITextFile(Phi.D2P(), BufferString, MPIOBJ);

  //===============================================================
  LogFile << "Finished Run \n" << std::endl;
  LogFile.close();
  MPI_Finalize();
  return 0;
}
