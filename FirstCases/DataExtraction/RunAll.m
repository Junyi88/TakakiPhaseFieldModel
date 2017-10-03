%% 0/9 Catalog

clear;
pname = 'E:\projects\TakakiPhaseFieldModel\FirstCases\G49 Strain10%'; %Location to extracte data and save data

NodesRaw.NodeFile=[pname '\nodes.inc'];
GetSetsFilename=[pname '\sets.inc'];
FilenameGetRotationMAtrix=[pname '\rotationMatrix.inc'];
SDVFilename=[pname '\SDVNodal.rpt'];
UFilename=[pname '\Displacements.rpt'];
numberOfNodes=211;
save names;

% NodesRaw.NodeFile='../G49/nodes.inc';

run ./AttemptToGetNodesAndEl.m;
run ./GetSets.m;
run ./GetRotationMAtrix.m;
run ./DealWithRptSDV.m;
run ./DealWithRptU.m;
run ./GLayer.m;
run ./DealWithRho.m;
run ./GetPhi.m;
run ./CombineData.m;