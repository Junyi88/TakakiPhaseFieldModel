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
#include "TakAngle.h"
#include "TakACBulkEnergy.h"
#include "TakACWallEnergy.h"
#include "TakACGradEnergy.h"
#include "TakACTOriEnergy.h"
#include "TakakiSolver.h"



// #include "GClass.h"

int main(int argc, char ** argv){

  const int NInputParameters=18;
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
    alpha, Wa, S, T, b, M0, Q, sg;

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
  T=InputParameters[11];
  InvPhiMin=InputParameters[12];
  b=InputParameters[13];
  M0=InputParameters[14];
  Q=InputParameters[15];
  sg = InputParameters[16];

  LogFile << "NY =  " << NY << std::endl;
  LogFile << "NX =  " << NX << std::endl;
  LogFile << "dy =  " << dy << std::endl; //delta y
  LogFile << "dx =  " << dx << std::endl; //delta x = 1 um
  LogFile << "dt =  " << dt << std::endl; 
  LogFile << "ntStart =  " << ntStart << std::endl; //start time
  LogFile << "ntEnd =  " << ntEnd << std::endl; //finish time
  LogFile << "WriteCount =  " << WriteCount << std::endl; //output interval
  LogFile << "MinAngle0 =  " << MinAngle0 << std::endl; //Delta theta_min
  LogFile << "BVec =  " << BVec << std::endl; //burgers vector
  LogFile << "mu =  " << mu << std::endl; //shear modulus MPa
  LogFile << "T =  " << T << std::endl; //temperature
  LogFile << "InvPhiMin =  " << InvPhiMin << std::endl;
  LogFile << "b =  " << b << std::endl; //
  LogFile << "M0 =  " << M0 << std::endl; //M_0
  LogFile << "Q =  " << Q << std::endl;  //
  LogFile << "sg =  " << sg << std::endl; //sigma
 

  // Calculate All
  const double kb = 1.38e-23;
  double M = M0*exp(-Q / (kb*T));
  
  Wa = 6.0*sg*b / (4.0*dx);
  alpha = sqrt(3.0*4.0*dx*sg / b);
  MTheta0 = M*sqrt(2 * Wa) / (6 * alpha);
  inMPhiConst = M*sqrt(2.0*Wa) / (6.0*alpha);
  S = sqrt(3.0*4.0*dx*sg / b);

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
  // InputParameters[11] 15=alpha
  // InputParameters[12] 16=Wa
  // InputParameters[13] 17=S

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
  Splitter(Eta0, BufferFull,  MPIOBJ);
  LogFile << "Reading Initial Conditions Eta Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Reading Initial Conditions Theta =  " << InitialConditionFileList[1] << std::endl;
  BufferString=InitialConditionFileList[1];
  ReadTextFile(BufferFull.Pointer() , BufferString, NX, NY, ',');
  Splitter(Theta0, BufferFull,  MPIOBJ);
  LogFile << "Reading Initial Conditions Theta Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Reading Initial Conditions Rho =  " << InitialConditionFileList[2] << std::endl;
  BufferString=InitialConditionFileList[2];
  ReadTextFile(BufferFull.Pointer() , BufferString, NX, NY, ',');
  Splitter( Rho0, BufferFull, MPIOBJ);
  LogFile << "Reading Initial Conditions Rho Completed \n" << std::endl;

  // *************************************************************************
  // FDClass=JFDMpi2DReflectHigh;
  // FDAngleClass=JFDMpi2DReflectHighAngle;
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

  LogFile << "========================================" << std::endl;
  LogFile << "Do 1 Step " << std::endl;
  Solver.Step_NoUpdate();
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

  LogFile << "Write Outputs - Theta" << std::endl;
  BufferString="Theta_0.csv";
  WriteMPITextFile(Theta.FP(), BufferString, MPIOBJ);

  BufferString="ThetaDx_0.csv";
  WriteMPITextFile(Theta.DxP(), BufferString, MPIOBJ);
  BufferString="ThetaDy_0.csv";
  WriteMPITextFile(Theta.DyP(), BufferString, MPIOBJ);
  BufferString="ThetaDxx_0.csv";
  WriteMPITextFile(Theta.DxxP(), BufferString, MPIOBJ);
  BufferString="ThetaDyy_0.csv";
  WriteMPITextFile(Theta.DyyP(), BufferString, MPIOBJ);
  BufferString="ThetaDxy_0.csv";
  WriteMPITextFile(Theta.DxyP(), BufferString, MPIOBJ);
  BufferString="ThetaD2_0.csv";
  WriteMPITextFile(Theta.D2P(), BufferString, MPIOBJ);

  LogFile << "Write Outputs - dFdPhase" << std::endl;
  BufferString="Bulk_dFdPhase_0.csv";
  WriteMPITextFile(BulkEnergy.dFdPhasePointer(), BufferString, MPIOBJ);
  BufferString="Grad_dFdPhase_0.csv";
  WriteMPITextFile(GradEnergy.dFdPhasePointer(), BufferString, MPIOBJ);
  BufferString="Wall_dFdPhase_0.csv";
  WriteMPITextFile(WallEnergy.dFdPhasePointer(), BufferString, MPIOBJ);
  BufferString="Ori_dFdPhase_0.csv";
  WriteMPITextFile(OriEnergy.dFdPhasePointer(), BufferString, MPIOBJ);

  LogFile << "Write Outputs - Others" << std::endl;
  BufferString="dThetadt_0.csv";
  WriteMPITextFile(OriEnergy.dThetadtPointer(), BufferString, MPIOBJ);
  BufferString="MTheta_0.csv";
  WriteMPITextFile(OriEnergy.MThetaPointer(), BufferString, MPIOBJ);

  BufferString="dPhidt_0.csv";
  WriteMPITextFile(Solver.dEtadtPointer(), BufferString, MPIOBJ);

  BufferString="Rho_0.csv";
  WriteMPITextFile(BulkEnergy.RhoPointer(), BufferString, MPIOBJ);
  BufferString="EStored_0.csv";
  WriteMPITextFile(BulkEnergy.EStoredPointer(), BufferString, MPIOBJ);

  BufferString=HeaderName + "_Phi_0.csv";
  WriteMPITextFile(Phi.FP(), BufferString, MPIOBJ);

  BufferString=HeaderName + "_Theta_0.csv";
  WriteMPITextFile(Theta.FP(), BufferString, MPIOBJ);
  // ***********************************************
  // Loop
  LogFile << "Start Loop \n" << std::endl;
  int counter=0;

  ntEnd++;

  for (int ntime=ntStart; ntime<ntEnd; ntime++){
    Solver.Step_All();
    if (counter<WriteCount){
      counter++;
    } else{
      counter=1;
      BufferString=HeaderName + "_Phi_" + std::to_string(ntime) + ".csv";
      WriteMPITextFile(Phi.FP(), BufferString, MPIOBJ);

      BufferString=HeaderName + "_Theta_" + std::to_string(ntime) + ".csv";
      WriteMPITextFile(Theta.FP(), BufferString, MPIOBJ);

      // BufferString=HeaderName + "_dEtadt_" + std::to_string(ntime) + ".csv";
      // WriteMPITextFile(Solver.dEtadtPointer(), BufferString, MPIOBJ);
      //
      // BufferString=HeaderName + "_dThetadt_" + std::to_string(ntime) + ".csv";
      // WriteMPITextFile(OriEnergy.dThetadtPointer(), BufferString, MPIOBJ);

      if (MPIOBJ.Nnode()==MPIOBJ.NLast())
        std::cout<<"nTime = "<< ntime <<std::endl;
    }

  }

  //===============================================================
  LogFile << "Finished Run \n" << std::endl;
  LogFile.close();
  MPI_Finalize();
  return 0;
}
