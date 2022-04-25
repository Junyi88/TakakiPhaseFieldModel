MOD_DIR=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SRC_DIR=$(MOD_DIR)src
INC_DIR=$(MOD_DIR)../inc
LIB_DIR=$(MOD_DIR)../obj
OBJ_DIR=$(MOD_DIR)../obj
EXEC_DIR=$(MOD_DIR)exec
CFLAGS=-std=c++11 -O2 -Wreorder

all : $(EXEC_DIR)/Test1.o $(EXEC_DIR)/Test2.o

$(EXEC_DIR)/Test1.o : $(OBJ_DIR)/JMat.o $(OBJ_DIR)/JMpi.o $(OBJ_DIR)/JFDMpi.o $(OBJ_DIR)/AuxFunctions.o $(SRC_DIR)/Test1.cpp
	mpicxx $(CFLAGS) -I$(INC_DIR) -I$(SRC_DIR) -o $(EXEC_DIR)/Test1.out $(SRC_DIR)/Test1.cpp \
		$(OBJ_DIR)/JMat.o $(OBJ_DIR)/JMpi.o $(OBJ_DIR)/JFDMpi.o $(OBJ_DIR)/AuxFunctions.o

$(EXEC_DIR)/Test2.o : $(OBJ_DIR)/JMat.o $(OBJ_DIR)/JMpi.o $(OBJ_DIR)/JFDMpi.o $(OBJ_DIR)/AuxFunctions.o $(SRC_DIR)/Test2.cpp
	mpicxx $(CFLAGS) -I$(INC_DIR) -I$(SRC_DIR) -o $(EXEC_DIR)/Test2.out $(SRC_DIR)/Test2.cpp \
		$(OBJ_DIR)/JMat.o $(OBJ_DIR)/JMpi.o $(OBJ_DIR)/JFDMpi.o $(OBJ_DIR)/AuxFunctions.o
