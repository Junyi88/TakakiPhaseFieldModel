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


figure(1);
clf;
hold on;
plot(X,'rx-');
plot(Y,'bs-');


[fitobject,gof,output] = fit([X,Y],Z,'linearinterp');

xMin=min(X);
xMax=max(X);
yMin=min(Y);
yMax=max(Y);

xx=linspace(xMin,xMax,50);
yy=linspace(yMin,yMax,50);
[XX,YY]=meshgrid(xx,yy);

ZZ=zeros(size(XX));

for ny=1:length(yy);
    for nx=1:length(xx);
        ZZ(ny,nx)=feval(fitobject,[XX(ny,nx),YY(ny,nx)]);
    end
end

figure(2);
clf;
hold on;
plot3(X,Y,Z,'rx');
surf(XX,YY,ZZ);