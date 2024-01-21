#include "AuxFunctions.h"

// ## ----------------------------------------------------------------
void ReadTextFile(double *ArrayPtr, std::string Filename, int NX, int NY, char delimiter)
{

  std::ifstream file;
	file.open(Filename);
	if (!file)
		std::cout << "Unable to open file: " << Filename <<std::endl;

    for (int i = 0; i < NY; ++i)
    {
        std::string line;
        std::getline(file, line);
        if ( !file.good() )
            break;

        std::stringstream iss(line);

        for (int j = 0; j < NX; ++j)
        {
            std::string val;
            std::getline(iss, val, delimiter);
            std::stringstream convertor(val);
            convertor >> *ArrayPtr;
            ArrayPtr++;
        }
    }

    file.close();
}

// ## ----------------------------------------------------------------
void Splitter(JMat & Part, JMat &Full, JMpi & MpiIn){

  for (int j=0; j<MpiIn.NYLo(); j++)
      for (int i=0; i<MpiIn.NX(); i++)
      {
        Part(j,i)=Full.Value(MpiIn.YStart()+j,i);
      }
}

// ## ----------------------------------------------------------------
void WriteMPITextFile(JMat & MatIn, const std::string &filename, JMpi & MpiIn){
  if (MpiIn.Nnode()==0) {
		std::ofstream file;
		file.open(filename, std::ofstream::out | std::ofstream::trunc);
    for (int j=0; j<MpiIn.NYLo(); j++){
      for (int i=0; i<(MpiIn.NX()-1); i++)
        file<<MatIn(j,i)<<',';
      file<<MatIn(j,(MpiIn.NX()-1))<<std::endl;

    }
    MPI_Send(MpiIn.WaitFlagPointer(),1,MPI_INT,MpiIn.Nnode()+1,81,MPI_COMM_WORLD);
    file.close();

  } else if (MpiIn.Nnode()==MpiIn.NLast()) {
    MPI_Recv(MpiIn.WaitFlagPointer(),1,MPI_INT,MpiIn.Nnode()-1,81,MPI_COMM_WORLD,MpiIn.MpiStatusPointer());
		std::ofstream file;
		file.open(filename, std::ofstream::out | std::ofstream::app);
    for (int j=0; j<MpiIn.NYLo(); j++){
      for (int i=0; i<(MpiIn.NX()-1); i++)
        file<<MatIn(j,i)<<',';
      file<<MatIn(j,(MpiIn.NX()-1))<<std::endl;
    }
		file.close();

  } else {
    MPI_Recv(MpiIn.WaitFlagPointer(),1,MPI_INT,MpiIn.Nnode()-1,81,MPI_COMM_WORLD,MpiIn.MpiStatusPointer());
		std::ofstream file;
		file.open(filename, std::ofstream::out | std::ofstream::app);
    for (int j=0; j<MpiIn.NYLo(); j++){
      for (int i=0; i<(MpiIn.NX()-1); i++)
        file<<MatIn(j,i)<<',';
      file<<MatIn(j,(MpiIn.NX()-1))<<std::endl;

    }
    MPI_Send(MpiIn.WaitFlagPointer(),1,MPI_INT,MpiIn.Nnode()+1,81,MPI_COMM_WORLD);
		file.close();
  }

}

// ## ----------------------------------------------------------------
void WriteMPITextFile(JMat * MatIn, const std::string &filename, JMpi & MpiIn){
  if (MpiIn.Nnode()==0) {
		std::ofstream file;
		file.open(filename, std::ofstream::out | std::ofstream::trunc);
    for (int j=0; j<MpiIn.NYLo(); j++){
      for (int i=0; i<(MpiIn.NX()-1); i++)
        file<<MatIn->Value(j,i)<<',';
      file<<MatIn->Value(j,(MpiIn.NX()-1))<<std::endl;

    }
    MPI_Send(MpiIn.WaitFlagPointer(),1,MPI_INT,MpiIn.Nnode()+1,81,MPI_COMM_WORLD);
    file.close();

  } else if (MpiIn.Nnode()==MpiIn.NLast()) {
    MPI_Recv(MpiIn.WaitFlagPointer(),1,MPI_INT,MpiIn.Nnode()-1,81,MPI_COMM_WORLD,MpiIn.MpiStatusPointer());
		std::ofstream file;
		file.open(filename, std::ofstream::out | std::ofstream::app);
    for (int j=0; j<MpiIn.NYLo(); j++){
      for (int i=0; i<(MpiIn.NX()-1); i++)
        file<<MatIn->Value(j,i)<<',';
      file<<MatIn->Value(j,(MpiIn.NX()-1))<<std::endl;
    }
		file.close();

  } else {
    MPI_Recv(MpiIn.WaitFlagPointer(),1,MPI_INT,MpiIn.Nnode()-1,81,MPI_COMM_WORLD,MpiIn.MpiStatusPointer());
		std::ofstream file;
		file.open(filename, std::ofstream::out | std::ofstream::app);
    for (int j=0; j<MpiIn.NYLo(); j++){
      for (int i=0; i<(MpiIn.NX()-1); i++)
        file<<MatIn->Value(j,i)<<',';
      file<<MatIn->Value(j,(MpiIn.NX()-1))<<std::endl;

    }
    MPI_Send(MpiIn.WaitFlagPointer(),1,MPI_INT,MpiIn.Nnode()+1,81,MPI_COMM_WORLD);
		file.close();
  }

}
