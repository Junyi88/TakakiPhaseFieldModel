MOD_DIR=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SRC_DIR=$(MOD_DIR)src
INC_DIR=$(MOD_DIR)inc
OBJ_DIR=$(MOD_DIR)obj
EXE_DIR=$(MOD_DIR)exec
FD_DIR=$(MOD_DIR)../JunyiFDMpiCpp/TwoD
FDSRC_DIR=$(MOD_DIR)../JunyiFDMpiCpp/TwoD/src
FDINC_DIR=$(MOD_DIR)../JunyiFDMpiCpp/TwoD/inc
FDOBJ_DIR=$(MOD_DIR)../JunyiFDMpiCpp/TwoD/obj
CFLAGS=-std=c++11 -O2 -Wreorder

all : $(OBJ_DIR)/TakPhase.o $(OBJ_DIR)/TakAngle.o \
	$(FDOBJ_DIR)/JMat.o $(FDOBJ_DIR)/JMpi.o $(FDOBJ_DIR)/JFDMpi.o $(FDOBJ_DIR)/AuxFunctions.o \
	$(FDOBJ_DIR)/JFDMpiAngle.o

$(OBJ_DIR)/TakAngle.o : $(SRC_DIR)/TakAngle.cpp $(INC_DIR)/TakAngle.h $(FDOBJ_DIR)/JFDMpiAngle.o
	mpicxx -c $(CFLAGS) -I$(FDINC_DIR) -I$(FDSRC_DIR) -I$(INC_DIR) -I$(SRC_DIR) -o \
		$(OBJ_DIR)/TakAngle.o $(SRC_DIR)/TakAngle.cpp \
		$(FDOBJ_DIR)/JMat.o $(FDOBJ_DIR)/JMpi.o $(FDOBJ_DIR)/JFDMpi.o $(FDOBJ_DIR)/JFDMpiAngle.o

$(OBJ_DIR)/TakPhase.o : $(SRC_DIR)/TakPhase.cpp $(INC_DIR)/TakPhase.h $(FDOBJ_DIR)/JFDMpi.o
	mpicxx -c $(CFLAGS) -I$(FDINC_DIR) -I$(FDSRC_DIR) -I$(INC_DIR) -I$(SRC_DIR) -o \
		$(OBJ_DIR)/TakPhase.o $(SRC_DIR)/TakPhase.cpp \
		$(FDOBJ_DIR)/JMat.o $(FDOBJ_DIR)/JMpi.o $(FDOBJ_DIR)/JFDMpi.o

$(FDOBJ_DIR)/JFDMpiAngle.o  : $(FDOBJ_DIR)/JFDMpi.o $(FDSRC_DIR)/JFDMpiAngle.cpp $(FDINC_DIR)/JFDMpiAngle.h
	mpicxx -c $(CFLAGS) -I$(FDINC_DIR) -I$(FDSRC_DIR) -o $(FDOBJ_DIR)/JFDMpiAngle.o $(FDSRC_DIR)/JFDMpiAngle.cpp \
		$(FDOBJ_DIR)/JMat.o $(FDOBJ_DIR)/JMpi.o $(FDOBJ_DIR)/JFDMpi.o

$(FDOBJ_DIR)/AuxFunctions.o  : $(FDOBJ_DIR)/JMat.o $(FDOBJ_DIR)/JMpi.o $(FDSRC_DIR)/AuxFunctions.cpp $(FDINC_DIR)/AuxFunctions.h
	mpicxx -c $(CFLAGS) -I$(FDINC_DIR) -I$(FDSRC_DIR) -o $(FDOBJ_DIR)/AuxFunctions.o $(FDSRC_DIR)/AuxFunctions.cpp \
		$(FDOBJ_DIR)/JMat.o $(FDOBJ_DIR)/JMpi.o

$(FDOBJ_DIR)/JFDMpi.o  : $(FDOBJ_DIR)/JMat.o $(FDOBJ_DIR)/JMpi.o $(FDSRC_DIR)/JFDMpi.cpp $(FDINC_DIR)/JFDMpi.h
	mpicxx -c $(CFLAGS) -I$(FDINC_DIR) -I$(FDSRC_DIR) -o $(FDOBJ_DIR)/JFDMpi.o $(FDSRC_DIR)/JFDMpi.cpp \
		$(FDOBJ_DIR)/JMat.o $(FDOBJ_DIR)/JMpi.o

$(FDOBJ_DIR)/JMpi.o  : $(FDSRC_DIR)/JMpi.cpp $(FDINC_DIR)/JMpi.h
	mpicxx -c $(CFLAGS) -I$(FDINC_DIR) -I$(FDSRC_DIR) -o $(FDOBJ_DIR)/JMpi.o $(FDSRC_DIR)/JMpi.cpp

$(FDOBJ_DIR)/JMat.o  : $(FDSRC_DIR)/JMat.cpp $(FDINC_DIR)/JMat.h
	mpicxx -c $(CFLAGS) -I$(FDINC_DIR) -I$(FDSRC_DIR) -o $(FDOBJ_DIR)/JMat.o $(FDSRC_DIR)/JMat.cpp
