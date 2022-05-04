MOD_DIR=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SRC_DIR=$(MOD_DIR)src
INC_DIR=$(MOD_DIR)inc
OBJ_DIR=$(MOD_DIR)obj
EXE_DIR=$(MOD_DIR)exec
TESTSRC_DIR=$(MOD_DIR)TestSrc
FD_DIR=$(MOD_DIR)../JunyiFDMpiCpp/TwoD
FDSRC_DIR=$(MOD_DIR)../JunyiFDMpiCpp/TwoD/src
FDINC_DIR=$(MOD_DIR)../JunyiFDMpiCpp/TwoD/inc
FDOBJ_DIR=$(MOD_DIR)../JunyiFDMpiCpp/TwoD/obj
CFLAGS=-std=c++11 -O2 -Wreorder

all : $(EXE_DIR)/Test1.out $(EXE_DIR)/Test1Mini.out $(EXE_DIR)/Test1Con.out

$(EXE_DIR)/Test1Con.out : $(TESTSRC_DIR)/Test1Con.cpp
	mpicxx $(CFLAGS) -I$(FDINC_DIR) -I$(FDSRC_DIR) -I$(INC_DIR) -I$(SRC_DIR) -o \
		$(EXE_DIR)/Test1Con.out $(TESTSRC_DIR)/Test1Con.cpp \
		$(FDSRC_DIR)/JMat.cpp $(FDSRC_DIR)/JMpi.cpp $(FDSRC_DIR)/JFDMpi.cpp \
		$(FDSRC_DIR)/AuxFunctions.cpp $(FDSRC_DIR)/JFDMpiAngle.cpp

$(EXE_DIR)/Test1Mini.out : $(TESTSRC_DIR)/Test1Mini.cpp
	mpicxx $(CFLAGS) -I$(FDINC_DIR) -I$(FDSRC_DIR) -I$(INC_DIR) -I$(SRC_DIR) -o \
		$(EXE_DIR)/Test1Mini.out $(TESTSRC_DIR)/Test1Mini.cpp \
		$(FDSRC_DIR)/JMat.cpp $(FDSRC_DIR)/JMpi.cpp $(FDSRC_DIR)/JFDMpi.cpp \
		$(FDSRC_DIR)/AuxFunctions.cpp $(FDSRC_DIR)/JFDMpiAngle.cpp

$(EXE_DIR)/Test1.out : $(TESTSRC_DIR)/Test1.cpp
	mpicxx $(CFLAGS) -I$(FDINC_DIR) -I$(FDSRC_DIR) -I$(INC_DIR) -I$(SRC_DIR) -o \
		$(EXE_DIR)/Test1.out $(TESTSRC_DIR)/Test1.cpp \
		$(FDSRC_DIR)/JMat.cpp $(FDSRC_DIR)/JMpi.cpp $(FDSRC_DIR)/JFDMpi.cpp \
		$(FDSRC_DIR)/AuxFunctions.cpp $(FDSRC_DIR)/JFDMpiAngle.cpp
