clear;

pname = 'E:\projects\TakakiPhaseFieldModel\FirstCases\G49 Strain10%';



NodesRaw.NodeFile=[pname '\nodes.inc'];
GetSetsFilename=[pname '\sets.inc'];
FilenameGetRotationMAtrix=[pname '\rotationMatrix.inc'];
SDVFilename=[pname '\SDVNodal.rpt'];
UFilename=[pname '\Displacements.rpt'];
numberOfNodes=211;
save names;

% NodesRaw.NodeFile='../G49/nodes.inc';
% GetSetsFilename='../G49/sets.inc';
% FilenameGetRotationMAtrix='../G49/rotationMatrix.inc';
% SDVFilename='../G49/SDVNodal.rpt';
% UFilename='../G49/Displacements.rpt';
% numberOfNodes=211;
% save names;

run ./AttemptToGetNodesAndEl.m;
run ./GetSets.m;
run ./GetRotationMAtrix.m;
run ./DealWithRptSDV.m;
run ./DealWithRptU.m;
run ./GLayer.m;
run ./DealWithRho.m;
run ./GetPhi.m;
run ./CombineData.m;