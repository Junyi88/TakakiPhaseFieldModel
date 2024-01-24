#include <ChenYunExec.h>


int main(int argc, char ** argv){
    
    //================================================================
    // # Initialise the MPI and system
    const int NInputParameters=17;
    int NPrs, Nnode;
    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &NPrs);
    MPI_Comm_rank(MPI_COMM_WORLD, &Nnode);

    // Load Files
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

    BufferInputStream1.close();
    for (int i=0; i<InitialConditionFileList.size(); i++)
        LogFile << "Input Files [" << i << "] : " << InitialConditionFileList[i] <<std::endl;
    LogFile << "Reading Initial Files List Completed \n" << std::endl;

    //====================================================================
    LogFile << "========================================" << std::endl;
    LogFile << "Reading Input Parameters " << Nnode << std::endl;
    LogFile << "Input Parameter File =  " << InputFile << std::endl;
    int NY, NX, ntStart, ntEnd, WriteCount;
    double dy, dx, dt, MinAngle0,
        epsilon, alpha, omega, beta, mu, tau_phi, tau_theta, w;

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

    epsilon = InputParameters[9];
    alpha = InputParameters[10];
    omega = InputParameters[11];
    beta = InputParameters[12];
    mu = InputParameters[13];
    tau_phi = InputParameters[14];
    tau_theta = InputParameters[15];
    w = InputParameters[16];

    double inv_tau_phi = 1.0 / tau_phi;

    //====================================================================
    // Setup MPI
    JMpi MPIOBJ(Nnode, NPrs, NY, NX, dy, dx);
    JMat BufferFull(NY, NX), Eta0(MPIOBJ.NYLo(), NX),
        Theta0(MPIOBJ.NYLo(), NX), Rho0(MPIOBJ.NYLo(), NX);

    // LOAD DATA
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

    TakPhase<JFDMpi2DReflectHigh> Phi(MPIOBJ, Eta0);
    TakAngle<JFDMpi2DReflectHighAngle> Theta(MPIOBJ, Theta0, MinAngle0);

    //=== Some Codes
    LogFile << "Start Loop \n" << std::endl;
    int counter=0;

    ntEnd++;
    for (int ntime=ntStart; ntime<ntEnd; ntime++){

        Phi.Calc_All();
        Theta.Calc_All();
        CYOriRHS_Term1<JFDMpi2DReflectHigh, JFDMpi2DReflectHighAngle> OriRHS_1(
            &Phi, &Theta, alpha, MinAngle0, MPIOBJ
        );
        CYOriRHS_Term2<JFDMpi2DReflectHigh, JFDMpi2DReflectHighAngle> OriRHS_2(
            &Phi, &Theta, omega, MPIOBJ
        );
        ChenYunOriLHS_Q<JFDMpi2DReflectHighAngle> OriLHS_Q(
            &Theta, beta, mu, omega, MPIOBJ
        );
        CYOri_DThetaDT<JFDMpi2DReflectHigh, JFDMpi2DReflectHighAngle> Ori_dTheta_dt(
            &OriLHS_Q, &OriRHS_1, &OriRHS_2, tau_theta, MPIOBJ
        );

        //====
        CYBulkRHS_Term1<JFDMpi2DReflectHigh> BulkRHS_1(&Phi, epsilon, MPIOBJ);
        CYBulkRHS_Term2<JFDMpi2DReflectHigh> BulkRHS_2(&Phi, epsilon, w, MPIOBJ);
        CYBulkRHS_Term3<JFDMpi2DReflectHigh, JFDMpi2DReflectHighAngle> BulkRHS_3(
            &Phi, &Theta, alpha, MPIOBJ);
        CYBulkRHS_Term4<JFDMpi2DReflectHigh, JFDMpi2DReflectHighAngle> BulkRHS_4(
            &Phi, &Theta, omega, MPIOBJ);

        CYBulk_DPhiDT dPhi_dt(
            &BulkRHS_1,
            &BulkRHS_2,
            &BulkRHS_3,
            &BulkRHS_4,
            tau_phi,
            MPIOBJ
        );

        for (int j=0; j < MPIOBJ.NYLo(); j++)
            for (int i=0; i < MPIOBJ.NX(); i++)
            {
                Phi.Update_Eta(dPhi_dt.val(j,i), dt, j, i);
                Theta.Update_Theta(Ori_dTheta_dt.val(j,u), dt, j, i);
            }

        if (counter<WriteCount){
      counter++;
    } else{
      counter=1;
      BufferString=HeaderName + "_Phi_" + std::to_string(ntime) + ".csv";
      WriteMPITextFile(Phi.FP(), BufferString, MPIOBJ);

      BufferString=HeaderName + "_Theta_" + std::to_string(ntime) + ".csv";
      WriteMPITextFile(Theta.FP(), BufferString, MPIOBJ);

      BufferString=HeaderName + "_dPhidt_" + std::to_string(ntime) + ".csv";
      WriteMPITextFile(Solver.dEtadtPointer(), BufferString, MPIOBJ);

      BufferString=HeaderName + "_dThetadt_" + std::to_string(ntime) + ".csv";
      WriteMPITextFile(OriEnergy.dThetadtPointer(), BufferString, MPIOBJ);

      if (MPIOBJ.Nnode()==MPIOBJ.NLast())
        std::cout<<"nTime = "<< ntime <<std::endl;
    }

    }

    if (MPIOBJ.Nnode()==MPIOBJ.NLast())
        std::cout<< "DONE" <<std::endl;

    //===============================================================
    LogFile << "Finished Run \n" << std::endl;
    LogFile.close();
    MPI_Finalize();
    return 0;
}
