#ifndef AUXFUNCTIONS_H
#define AUXFUNCTIONS_H

#include <iostream>
#include <string>
#include <istream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <mpi.h>
#include "JMpi.h"
#include "JMat.h"

void ReadTextFile(double *ArrayPtr, std::string Filename, int NX, int NY, char delimiter);
void Splitter(JMat & Part, JMat &Full, JMpi & MpiIn);
void WriteMPITextFile(JMat & MatIn, const std::string &filename, JMpi & MpiIn);
void WriteMPITextFile(JMat * MatIn, const std::string &filename, JMpi & MpiIn);
#endif
