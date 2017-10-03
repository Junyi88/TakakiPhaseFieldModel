%% 6/9 to run: Get Displacement.rpt (i.e.U1 U2 U3_Displacemnt)

clear;
load NodeDataIgnore.mat;
load SDVDataIgnore.mat;
load DisplacementsDataIgnore.mat;
for n1=1:length(NodeLayers)
   SDV(n1).S=zeros(length(Nodes(n1).N),3); 
   SDV(n1).U=zeros(length(Nodes(n1).N),3); 
end

for n2=1:length(NodeLayers)
for n1=1:length((Nodes(n2).N))
    SDV(n2).S(n1,1:2)=Nodes(n2).N(n1,2:3);
    SDV(n2).S(n1,3)=SDVs(Nodes(n2).N(n1,1),4); %DANGER HARD CODING
    SDV(n2).U(n1,:)=Us(Nodes(n2).N(n1,1),2:4);
end
end

for n1=1:length(NodeLayers)
   SDV(n1).Xo=SDV(n1).S(:,1); 
   SDV(n1).Yo=SDV(n1).S(:,2); 
   SDV(n1).Zo=SDV(n1).S(:,3); 
   
   SDV(n1).X=SDV(n1).Xo;%+SDV(n2).U(n1,1); 
   SDV(n1).Y=SDV(n1).Yo;%+SDV(n2).U(n1,2); 
   SDV(n1).Z=SDV(n1).Zo; 
end

for n1=1:length(NodeLayers)
[SDV(n1).fitobject,SDV(n1).gof,SDV(n1).output] = fit([SDV(n1).X,SDV(n1).Y],SDV(n1).Z,'linearinterp');
end

xMin=min(SDV(1).Xo);
xMax=max(SDV(1).Xo);
yMin=min(SDV(1).Yo);
yMax=max(SDV(1).Yo);

xx=linspace(xMin,xMax,201);
yy=linspace(yMin,yMax,201);
[XX,YY]=meshgrid(xx,yy);

for n1=1:length(NodeLayers)
    SDV(n1).ZZ=zeros(size(XX));
for ny=1:length(yy);
    for nx=1:length(xx);
        SDV(n1).ZZ(ny,nx)=feval(SDV(n1).fitobject,[XX(ny,nx),YY(ny,nx)]);
        
    end
end
end

n1=20;
% figure(61);
% clf;
% hold on;
% plot(XX(n1,:),SDV(1).ZZ(n1,:),'rx-');
% plot(XX(n1,:),SDV(2).ZZ(n1,:),'gs-');
% plot(XX(n1,:),SDV(3).ZZ(n1,:),'b>-');
% plot(XX(n1,:),SDV(4).ZZ(n1,:),'ks-');
% plot(XX(n1,:),SDV(5).ZZ(n1,:),'cx-');
% plot(XX(n1,:),SDV(6).ZZ(n1,:),'y>-');
% plot(XX(n1,:),SDV(7).ZZ(n1,:),'mx-');
% plot(XX(n1,:),SDV(8).ZZ(n1,:),'rs-');
% plot(XX(n1,:),SDV(9).ZZ(n1,:),'cs-');
% 
% figure(62);
% clf;
% hold on;
% plot3(XX(:),YY(:),SDV(1).ZZ(:),'rx');
% plot3(XX(:),YY(:),SDV(2).ZZ(:),'gs');
% plot3(XX(:),YY(:),SDV(3).ZZ(:),'b>');
% plot3(XX(:),YY(:),SDV(4).ZZ(:),'ks');
% plot3(XX(:),YY(:),SDV(5).ZZ(:),'cx');
% plot3(XX(:),YY(:),SDV(6).ZZ(:),'y>');
% plot3(XX(:),YY(:),SDV(7).ZZ(:),'mx');
% plot3(XX(:),YY(:),SDV(8).ZZ(:),'rs');
% plot3(XX(:),YY(:),SDV(9).ZZ(:),'cs');

save SDVSortedIgnore SDV;