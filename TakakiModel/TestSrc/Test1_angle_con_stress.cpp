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
#include "BasicChemPotential.h"

#include "TakACBulkEnergy.h"
#include "TakACWallEnergy.h"
#include "TakACGradEnergy.h"
#include "TakACTOriEnergyCon.h"

#include "TakACChemEnergy.h"

#include "TakakiSolverAngleCon.h"

int main(int argc, char ** argv){

  const int NInputParameters=23;
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

  BufferInputStream1 >> BufferString;
  InitialConditionFileList.push_back(BufferString); // Con

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
    alpha, Wa, S, T, b, M0, Q, sg, 
    kappa_chem, MChem;

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

  kappa_chem = InputParameters[17];
  MChem = InputParameters[18];

  double Stress = InputParameters[19];
  double MPhiStress = InputParameters[20];
  double MThetaStress = InputParameters[21];
  double MChemStress = InputParameters[22];


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

  LogFile << "kappa_chem =  " << kappa_chem << std::endl; //sigma
  LogFile << "MChem =  " << MChem << std::endl; //sigma

  LogFile << "Stress =  " << Stress << std::endl; //sigma
  LogFile << "MPhiStress =  " << MPhiStress << std::endl; //sigma
  LogFile << "MThetaStress =  " << MThetaStress << std::endl; //sigma
  LogFile << "MChemStress =  " << MChemStress << std::endl; //sigma

  if (Stress > 0.0)
  {
    if (MPhiStress > 0.0){
      inMPhiConst *= Stress * MPhiStress;
    }
    if (MThetaStress > 0.0){
      MTheta0 *= Stress * MThetaStress;
    }
    if (MChemStress > 0.0){
      MChem *= Stress * MChemStress;
    }
  }
    


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

  JMat Con0(MPIOBJ.NYLo(), NX);

  LogFile << "Setup MPI Object Completed \n" << std::endl;

  LogFile << "Nnode: " << MPIOBJ.Nnode() << std::endl;
  LogFile << "NPrs: " << MPIOBJ.NPrs() << std::endl;
  LogFile << "Nwork: " << MPIOBJ.Nwork() << std::endl;
  LogFile << "NLast: " << MPIOBJ.NLast() << std::endl;
  LogFile << "dy: " << MPIOBJ.dy() << std::endl;
  LogFile << "dx: " << MPIOBJ.dx() << std::endl;
  LogFile << "NX: " << MPIOBJ.NX() << std::endl;

  LogFile << "NYFu: " << MPIOBJ.NYFu() << std::endl;
  LogFile << "NYLo: " << MPIOBJ.NYLo() << std::endl;
  LogFile << "NYGl: " << MPIOBJ.NYGl() << std::endl;
  LogFile << "NYShort: " << MPIOBJ.NYShort() << std::endl;
  LogFile << "NXYGl: " << MPIOBJ.NXYGl() << std::endl;
  LogFile << "PosShortLast: " << MPIOBJ.PosShortLast() << std::endl;
  LogFile << "YStart: " << MPIOBJ.YStart() << std::endl;
  // LogFile << "Nnode: " << JMpi.Nnode() << std::endl;



  LogFile << "========================================" << std::endl;
  LogFile << "Reading Initial Conditions Eta =  " << InitialConditionFileList[0] << std::endl;
  BufferString=InitialConditionFileList[0];
  ReadTextFile(BufferFull.Pointer() , BufferString, NX, NY, ',');
  LogFile << "Splitter \n" << std::endl;
  Splitter(Eta0, BufferFull, MPIOBJ);
  LogFile << std::endl;
  LogFile << "Reading Initial Conditions Eta Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Reading Initial Conditions Theta =  " << InitialConditionFileList[1] << std::endl;
  BufferString=InitialConditionFileList[1];
  ReadTextFile(BufferFull.Pointer() , BufferString, NX, NY, ',');
  Splitter(Theta0, BufferFull, MPIOBJ);
  LogFile << "Reading Initial Conditions Theta Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Reading Initial Conditions Rho =  " << InitialConditionFileList[2] << std::endl;
  BufferString=InitialConditionFileList[2];
  ReadTextFile(BufferFull.Pointer() , BufferString, NX, NY, ',');
  Splitter(Rho0, BufferFull, MPIOBJ);
  LogFile << "Reading Initial Conditions Rho Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Reading Initial Conditions Con =  " << InitialConditionFileList[3] << std::endl;
  BufferString=InitialConditionFileList[3];
  ReadTextFile(BufferFull.Pointer() , BufferString, NX, NY, ',');
  Splitter(Con0, BufferFull, MPIOBJ);
  LogFile << "Reading Initial Conditions Con Completed \n" << std::endl;

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
  LogFile << "Setup Con Class " << std::endl;
  BasicChemPotential<JFDMpi2DExternal1> Con(MPIOBJ, kappa_chem, Con0);
  LogFile << "Setup Con Class Completed \n" << std::endl;

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
  TakACTOriEnergyCon<JFDMpi2DReflectHigh, JFDMpi2DReflectHighAngle, JFDMpi2DExternal1> OriEnergy(&Phi, &Theta, &Con,
    S, MTheta0, InvPhiMin, MPIOBJ);
  LogFile << "Setup TakACTOriEnergy Class Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Setup TakACChemEnergy Class " << std::endl;
  TakACChemEnergy<JFDMpi2DExternal1> ChemEnegy(&Con, MChem, MPIOBJ);
  LogFile << "Setup TakACChemEnergy Class Completed \n" << std::endl;

  LogFile << "========================================" << std::endl;
  LogFile << "Setup TakakiSolver Class " << std::endl;
  TakakiSolverAngleCon<JFDMpi2DReflectHigh, JFDMpi2DReflectHighAngle, JFDMpi2DExternal1> Solver(
    &Phi, &Theta, &Con,
    &BulkEnergy, &WallEnergy, &GradEnergy, &OriEnergy, &ChemEnegy,
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

  BufferString="Phi_0_R.csv";
  WriteMPITextFile(Phi.FP(), BufferString, MPIOBJ);

  BufferString=HeaderName + "_Con_0_R.csv";
  WriteMPITextFile(Con.FP(), BufferString, MPIOBJ);

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

  BufferString=HeaderName + "_Con_0.csv";
  WriteMPITextFile(Con.FP(), BufferString, MPIOBJ);

//  BufferString=HeaderName + "_Q1_0.csv";
//  WriteMPITextFile(Theta.FP1(), BufferString, MPIOBJ);
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


      BufferString=HeaderName + "_Con_" + std::to_string(ntime) + ".csv";
      WriteMPITextFile(Con.FP(), BufferString, MPIOBJ);
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
