#include <iostream>
#include <string>
#include <istream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>

#include "JMpi.h"
#include "JMat.h"
#include "JFDMpi.h"
#include "JFDMpiAngle.h"
#include "AuxFunctions.h"

#include "TakPhase.h"
#include "TakAngle.h"
#include "TakACBulkEnergy.h"
#include "TakACWallEnergy.h"
#include "TakACGradEnergy.h"
#include "TakACTOriEnergy.h"
#include "TakakiSolver.h"

int main(int argc, char ** argv){

  const int NInputParameters=17;
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

  InitialConditionFile=argv[2];
  InputFile=argv[3];

  //====================================================================
  LogFile << "========================================" << std::endl;
  LogFile << "Reading Initial Files List On Node" << Nnode << std::endl;
  BufferInputStream1.open(InitialConditionFile);

  BufferInputStream1>>BufferString;
  InitialConditionFileList.push_back(BufferString); // Eta
  BufferInputStream1>>BufferString;
  InitialConditionFileList.push_back(BufferString); // Theta
  BufferInputStream1>>BufferString;
  InitialConditionFileList.push_back(BufferString); // Rho

  BufferInputStream1.close();
  for (int i=0; i<InitialConditionFileList.size(); i++)
    LogFile << "Input Files [" << i << "] : " << InitialConditionFileList[i] <<std::endl;
  LogFile << "Reading Initial Files List Completed \n" << std::endl;

  //====================================================================
  LogFile << "========================================" << std::endl;
  LogFile << "Reading Input Parameters " << Nnode << std::endl;
  LogFile << "Input Parameter File =  " << InputFile << std::endl;
  int NY, NX, ntStart, ntEnd, WriteCount;
  double dy, dx, dt, MinAngle0, BVec, mu, MTheta0, InvPhiMin, inMPhiConst,
    alpha, Wa, S;

  double InputParameters[50];
  ReadTextFile(&InputParameters[0], InputFile, 1, NInputParameters, ',');

  NY=InputParameters[0];
  NX=InputParameters[1];
  dy=InputParameters[2];
  dx=InputParameters[3];
  dt=InputParameters[4];
  ntStart=InputParameters[5];
  ntEnd=InputParameters[6];
  WriteCount=InputParameters[7];
  MinAngle0=InputParameters[8];
  BVec=InputParameters[9];
  mu=InputParameters[10];
  MTheta0=InputParameters[11];
  InvPhiMin=InputParameters[12];
  inMPhiConst=InputParameters[13];
  alpha=InputParameters[14];
  Wa=InputParameters[15];
  S=InputParameters[16];

  LogFile << "NY =  " << NY << std::endl;
  LogFile << "NX =  " << NX << std::endl;
  LogFile << "dy =  " << dy << std::endl;
  LogFile << "dx =  " << dx << std::endl;
  LogFile << "dt =  " << dt << std::endl;
  LogFile << "ntStart =  " << ntStart << std::endl;
  LogFile << "ntEnd =  " << ntEnd << std::endl;
  LogFile << "WriteCount =  " << WriteCount << std::endl;
  LogFile << "MinAngle0 =  " << MinAngle0 << std::endl;
  LogFile << "BVec =  " << BVec << std::endl;
  LogFile << "mu =  " << mu << std::endl;
  LogFile << "MTheta0 =  " << MTheta0 << std::endl;
  LogFile << "InvPhiMin =  " << InvPhiMin << std::endl;
  LogFile << "inMPhiConst =  " << inMPhiConst << std::endl;
  LogFile << "alpha =  " << alpha << std::endl;
  LogFile << "Wa =  " << Wa << std::endl;
  LogFile << "S =  " << S << std::endl;
  // Input Parameters List
  // InputParameters[0] 1=NY
  // InputParameters[1] 2=NX
  // InputParameters[2] 3=dy
  // InputParameters[3] 4=dx
  // InputParameters[4] 5=dt
  // InputParameters[5] 6=ntStart
  // InputParameters[6] 7=ntEnd
  // InputParameters[7] 8=WriteCount
  // InputParameters[8] 9=MinAngle0
  // InputParameters[9] 10=BVec
  // InputParameters[10] 11=mu
  // InputParameters[11] 12=MTheta0
  // InputParameters[12] 13=InvPhiMin
  // InputParameters[13] 14=inMPhiConst

  LogFile << "Reading Parameter Files Completed \n" << std::endl;

  //====================================================================
  LogFile << "========================================" << std::endl;
  LogFile << "Setup MPI Object" << std::endl;

  JMpi MPIOBJ(Nnode, NPrs, NY, NX, dy, dx);
  JMat BufferFull(NY, NX), Eta0(MPIOBJ.NYLo(), NX),
    Theta0(MPIOBJ.NYLo(), NX), Rho0(MPIOBJ.NYLo(), NX);

  LogFile << "Setup MPI Object Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Reading Initial Conditions Eta =  " << InitialConditionFileList[0] << std::endl;
  BufferString=InitialConditionFileList[0];
  ReadTextFile(BufferFull.Pointer() , BufferString, NX, NY, ',');
  Splitter(BufferFull, Eta0, MPIOBJ);
  LogFile << "Reading Initial Conditions Eta Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Reading Initial Conditions Theta =  " << InitialConditionFileList[1] << std::endl;
  BufferString=InitialConditionFileList[1];
  ReadTextFile(BufferFull.Pointer() , BufferString, NX, NY, ',');
  Splitter(BufferFull, Theta0, MPIOBJ);
  LogFile << "Reading Initial Conditions Theta Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Reading Initial Conditions Rho =  " << InitialConditionFileList[2] << std::endl;
  BufferString=InitialConditionFileList[2];
  ReadTextFile(BufferFull.Pointer() , BufferString, NX, NY, ',');
  Splitter(BufferFull, Rho0, MPIOBJ);
  LogFile << "Reading Initial Conditions Rho Completed \n" << std::endl;

  //*************************************************************************
  //FDClass=JFDMpi2DReflectHigh;
  //FDAngleClass=JFDMpi2DReflectHighAngle;
  LogFile << "========================================" << std::endl;
  LogFile << "Setup Phase Class " << std::endl;
  TakPhase<JFDMpi2DReflectHigh> Phi(MPIOBJ, Eta0);
  LogFile << "Setup Phase Class Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Setup Theta Class " << std::endl;
  TakAngle<JFDMpi2DReflectHighAngle> Theta(MPIOBJ, Theta0, MinAngle0);
  LogFile << "Setup Theta Class Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Setup TakACBulkEnergy Class " << std::endl;
  TakACBulkEnergy<JFDMpi2DReflectHigh> BulkEnergy(&Phi, BVec, mu, Rho0, MPIOBJ);
  LogFile << "Setup TakACBulkEnergy Class Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Setup TakACGradEnergy Class " << std::endl;
  TakACGradEnergy<JFDMpi2DReflectHigh> GradEnergy(&Phi, alpha, MPIOBJ);
  LogFile << "Setup TakACGradEnergy Class Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Setup TakACGradEnergy Class " << std::endl;
  TakACWallEnergy<JFDMpi2DReflectHigh> WallEnergy(&Phi, Wa, MPIOBJ);
  LogFile << "Setup TakACGradEnergy Class Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Setup TakACTOriEnergy Class " << std::endl;
  TakACTOriEnergy<JFDMpi2DReflectHigh, JFDMpi2DReflectHighAngle> OriEnergy(&Phi, &Theta,
    S, MTheta0, InvPhiMin, MPIOBJ);
  LogFile << "Setup TakACTOriEnergy Class Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Setup TakakiSolver Class " << std::endl;
  TakakiSolver<JFDMpi2DReflectHigh, JFDMpi2DReflectHighAngle> Solver(&Phi, &Theta,
    &BulkEnergy, &WallEnergy, &GradEnergy, &OriEnergy,
    inMPhiConst , dt, MPIOBJ);
  LogFile << "Setup TakakiSolver Class Completed \n" << std::endl;
  //TakPhase(JMpi inJMpi, const JMat &in_F);
  //TakAngle(JMpi inJMpi, const JMat &in_F, const double &inMinMag);
  // TakACBulkEnergy(TakPhase<FDClass> * inPhi, const double &BVec, const double &mu,
  //   const JMat &RhoIn, JMpi inJMpi);
  // TakACGradEnergy(TakPhase<FDClass> * inPhi, const double &inalpha, JMpi inJMpi);
  // TakACWallEnergy(TakPhase<FDClass> * inPhi, const double &inWa, JMpi inJMpi);
  // TakACTOriEnergy(TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta,
  //   const double &insConst, const double &inMTheta0,  const double &inInvPhiMin,
  //   JMpi inJMpi);
  // TakakiSolver(TakPhase<FDClass> * inPhi, TakAngle<FDAngleClass> * inTheta,
  //   TakACBulkEnergy<FDClass> * inBulkEnergy, TakACWallEnergy<FDClass> * inWallEnergy,
  //   TakACGradEnergy<FDClass> * inGradEnergy, TakACTOriEnergy<FDClass, FDAngleClass> * inOriEnergy,
  //   const double &inMPhiConst, const double &indt,
  //   JMpi inJMpi);

  LogFile << "========================================" << std::endl;
  LogFile << "Do 1 Step " << std::endl;
  Solver.Step_NoUpdate();
  LogFile << "Done \n" << std::endl;

  //===============================================================
  LogFile << "Finished Run \n" << std::endl;
  LogFile.close();
  MPI_Finalize();
  return 0;
}
