clear;
load NodeDataIgnore.mat;
load SDVDataIgnore.mat;


SDVLayer1=zeros(length(Nodes(1).N),3);
for n1=1:length((Nodes(1).N))
    SDVLayer1(n1,1:2)=Nodes(1).N(n1,2:3);
    SDVLayer1(n1,3)=SDVs(Nodes(1).N(n1,1),3);
end

X=SDVLayer1(:,1);
Y=SDVLayer1(:,2);
Z=SDVLayer1(:,3);


SDVLayer2=zeros(length(Nodes(5).N),3);
for n1=1:length((Nodes(5).N))
    SDVLayer2(n1,1:2)=Nodes(5).N(n1,2:3);
    SDVLayer2(n1,3)=SDVs(Nodes(5).N(n1,1),3);
end

X2=SDVLayer2(:,1);
Y2=SDVLayer2(:,2);
Z2=SDVLayer2(:,3);

figure(1);
clf;
hold on;
plot(X,'rx-');
plot(Y,'bs-');


[fitobject,gof,output] = fit([X,Y],Z,'linearinterp');
[fitobject2,gof2,output2] = fit([X2,Y2],Z2,'linearinterp');
xMin=min(X);
xMax=max(X);
yMin=min(Y);
yMax=max(Y);

xx=linspace(xMin,xMax,50);
yy=linspace(yMin,yMax,50);
[XX,YY]=meshgrid(xx,yy);

ZZ=zeros(size(XX));
ZZ2=zeros(size(XX));
for ny=1:length(yy);
    for nx=1:length(xx);
        ZZ(ny,nx)=feval(fitobject,[XX(ny,nx),YY(ny,nx)]);
        ZZ2(ny,nx)=feval(fitobject2,[XX(ny,nx),YY(ny,nx)]);
    end
end

figure(2);
clf;
hold on;
% plot3(X,Y,Z,'rx');
surf(XX,YY,ZZ);
plot3(XX(:),YY(:),ZZ2(:),'rx');

figure(3);
clf;
hold on;
% plot3(X,Y,Z,'rx');
surf(XX,YY,abs(ZZ-ZZ2)./abs(ZZ));
% plot3(XX(:),YY(:),ZZ2(:),'rx');