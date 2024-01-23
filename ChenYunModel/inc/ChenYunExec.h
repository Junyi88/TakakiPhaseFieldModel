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

#include <ChenYunOriRHS_Term1.h>
#include <ChenYunOriRHS_Term2.h>
#include <ChenYunOriLHS_Q.h>
#include <ChenYunDThetaDT.h>

#include <ChenYunBulkRHS_Term1.h>
#include <ChenYunBulkRHS_Term2.h>
#include <ChenYunBulkRHS_Term3.h>
#include <ChenYunBulkRHS_Term4.h>
#include <CYBulk_DPhiDT.h>
