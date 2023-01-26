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
  double MChem = 1.0;
  double kappacon = 0.01;

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
  Splitter(Con0, BufferFull,  MPIOBJ);
  LogFile << "Reading Initial Conditions Eta Completed \n" << std::endl;

  // *********************************X****************************************
  // FDClass=JFDMpi2DReflectHigh;
  // FDAngleClass=JFDMpi2DReflectHighAngle;
  LogFile << "========================================" << std::endl;
  LogFile << "Setup Con Class " << std::endl;
  BasicChemPotential<JFDMpi2DExternal1> Con(MPIOBJ, kappacon, Con0);
  LogFile << "Setup Con Class Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Setup TakACChemEnergy Class " << std::endl;
  TakACChemEnergy<JFDMpi2DExternal1> ChemEnegy(&Con, MChem, MPIOBJ);
  LogFile << "Setup TakACChemEnergy Class Completed \n" << std::endl;


  // LogFile << "========================================" << std::endl;
  // LogFile << "Setup TakACChemEnergy Class " << std::endl;
  // TakACChemEnergy<JFDMpi2DReflectHigh> ChemEnegy(&Con, MChem, MPIOBJ);
  // LogFile << "Setup TakACChemEnergy Class Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Do 1 Step " << std::endl;
  Con.Calc_All();
  ChemEnegy.Calc_All();
  LogFile << "Done \n" << std::endl;

  //===============================================================
  // Write outputs For First STep
  LogFile << "========================================" << std::endl;
  LogFile << "Write Outputs - Con" << std::endl;

  BufferString="Con_0.csv";
  WriteMPITextFile(Con.FP(), BufferString, MPIOBJ);
  BufferString="Mu_0.csv";
  WriteMPITextFile(Con.Mu(), BufferString, MPIOBJ);

  BufferString="D2Mu_0.csv";
  WriteMPITextFile(Con.DMuPtr()->D2P(), BufferString, MPIOBJ);

  BufferString="dcondt_0.csv";
  WriteMPITextFile(ChemEnegy.dcondt(), BufferString, MPIOBJ);

  double _Ny = MPIOBJ.NYLo();
  double _NX = MPIOBJ.NX();
  double dtime = 0.001;

  for (int j = 0; j<_Ny; j++)
    for (int i = 0; i<_NX; i++) {
      Con.Update_Con(ChemEnegy.dcondt()->Value(j, i), dtime, j, i);
    }
  BufferString="Con_1.csv";
  WriteMPITextFile(Con.FP(), BufferString, MPIOBJ);

  for (int ntime = 2; ntime < 2000; ntime++)
  {
      Con.Calc_All();
      ChemEnegy.Calc_All();

      for (int j = 0; j<_Ny; j++)
        for (int i = 0; i<_NX; i++) {
          Con.Update_Con(ChemEnegy.dcondt()->Value(j, i), dtime, j, i);
        }

      BufferString="Con_" + std::to_string(ntime) + ".csv";
      WriteMPITextFile(Con.FP(), BufferString, MPIOBJ);
      BufferString="Mu_0" + std::to_string(ntime) + ".csv";
      WriteMPITextFile(Con.Mu(), BufferString, MPIOBJ);

      BufferString="D2Mu_" + std::to_string(ntime) + ".csv";
      WriteMPITextFile(Con.DMuPtr()->D2P(), BufferString, MPIOBJ);

      BufferString="dcondt_" + std::to_string(ntime) + ".csv";
      WriteMPITextFile(ChemEnegy.dcondt(), BufferString, MPIOBJ);


  }

  // BufferString=HeaderName + "_Theta_" + std::to_string(ntime) + ".csv";

  // BufferString="Phi3_0.csv";
  // WriteMPITextFile(Phi.F3P(), BufferString, MPIOBJ);
  // BufferString="Phi4_0.csv";
  // WriteMPITextFile(Phi.F4P(), BufferString, MPIOBJ);
  // BufferString="Phi5_0.csv";
  // WriteMPITextFile(Phi.F5P(), BufferString, MPIOBJ);

  // BufferString="PhiDx_0.csv";
  // WriteMPITextFile(Phi.DxP(), BufferString, MPIOBJ);
  // BufferString="PhiDy_0.csv";
  // WriteMPITextFile(Phi.DyP(), BufferString, MPIOBJ);
  // BufferString="PhiDxx_0.csv";
  // WriteMPITextFile(Phi.DxxP(), BufferString, MPIOBJ);
  // BufferString="PhiDyy_0.csv";
  // WriteMPITextFile(Phi.DyyP(), BufferString, MPIOBJ);
  // BufferString="PhiDxy_0.csv";
  // WriteMPITextFile(Phi.DxyP(), BufferString, MPIOBJ);
  // BufferString="PhiD2_0.csv";
  // WriteMPITextFile(Phi.D2P(), BufferString, MPIOBJ);

  //===============================================================
  LogFile << "Finished Run \n" << std::endl;
  LogFile.close();
  MPI_Finalize();
  return 0;
}
