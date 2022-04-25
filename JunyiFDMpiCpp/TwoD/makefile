MOD_DIR=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SRC_DIR=$(MOD_DIR)src
INC_DIR=$(MOD_DIR)inc
LIB_DIR=$(MOD_DIR)lib
OBJ_DIR=$(MOD_DIR)obj
CFLAGS=-std=c++11 -O2 -Wreorder

all : $(OBJ_DIR)/JMat.o $(OBJ_DIR)/JMpi.o $(OBJ_DIR)/JFDMpi.o $(OBJ_DIR)/AuxFunctions.o \
	$(OBJ_DIR)/JFDMpiAngle.o

$(OBJ_DIR)/JFDMpiAngle.o  : $(OBJ_DIR)/JFDMpi.o $(SRC_DIR)/JFDMpiAngle.cpp $(INC_DIR)/JFDMpiAngle.h
	mpicxx -c $(CFLAGS) -I$(INC_DIR) -I$(SRC_DIR) -o $(OBJ_DIR)/JFDMpiAngle.o $(SRC_DIR)/JFDMpiAngle.cpp \
		$(OBJ_DIR)/JMat.o $(OBJ_DIR)/JMpi.o $(OBJ_DIR)/JFDMpi.o

$(OBJ_DIR)/AuxFunctions.o  : $(OBJ_DIR)/JMat.o $(OBJ_DIR)/JMpi.o $(SRC_DIR)/AuxFunctions.cpp $(INC_DIR)/AuxFunctions.h
	mpicxx -c $(CFLAGS) -I$(INC_DIR) -I$(SRC_DIR) -o $(OBJ_DIR)/AuxFunctions.o $(SRC_DIR)/AuxFunctions.cpp \
		$(OBJ_DIR)/JMat.o $(OBJ_DIR)/JMpi.o

$(OBJ_DIR)/JFDMpi.o  : $(OBJ_DIR)/JMat.o $(OBJ_DIR)/JMpi.o $(SRC_DIR)/JFDMpi.cpp $(INC_DIR)/JFDMpi.h
	mpicxx -c $(CFLAGS) -I$(INC_DIR) -I$(SRC_DIR) -o $(OBJ_DIR)/JFDMpi.o $(SRC_DIR)/JFDMpi.cpp \
		$(OBJ_DIR)/JMat.o $(OBJ_DIR)/JMpi.o

$(OBJ_DIR)/JMpi.o  : $(SRC_DIR)/JMpi.cpp $(INC_DIR)/JMpi.h
	mpicxx -c $(CFLAGS) -I$(INC_DIR) -I$(SRC_DIR) -o $(OBJ_DIR)/JMpi.o $(SRC_DIR)/JMpi.cpp

$(OBJ_DIR)/JMat.o  : $(SRC_DIR)/JMat.cpp $(INC_DIR)/JMat.h
	mpicxx -c $(CFLAGS) -I$(INC_DIR) -I$(SRC_DIR) -o $(OBJ_DIR)/JMat.o $(SRC_DIR)/JMat.cpp
