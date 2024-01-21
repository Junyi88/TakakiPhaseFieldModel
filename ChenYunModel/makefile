MOD_DIR=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SRC_DIR=$(MOD_DIR)src
INC_DIR=$(MOD_DIR)inc
OBJ_DIR=$(MOD_DIR)obj
EXE_DIR=$(MOD_DIR)exec

CFLAGS=-std=c++11 -O2 -Wreorder -march=native

all : $(EXE_DIR)/ChenYunExec.out

$(EXE_DIR)/ChenYunExec.out : $(SRC_DIR)/ChenYunExec.cpp
	mpicxx $(CFLAGS) -I$(INC_DIR) -o \
		$(EXE_DIR)/ChenYunExec.out $(SRC_DIR)/ChenYunExec.cpp \
		$(SRC_DIR)/JMat.cpp $(SRC_DIR)/JMpi.cpp $(SRC_DIR)/JFDMpi.cpp \
		$(SRC_DIR)/AuxFunctions.cpp $(SRC_DIR)/JFDMpiAngle.cpp
		

.PHONY: clean
clean:
		rm -f $(EXE_DIR)/*.out
