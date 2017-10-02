clear;

load SetsIgnore;
load EulerAnglesIgnore;
load SDVSortedIgnore;
load NodeDataIgnore;

%%
load names;

Nnodes=numberOfNodes;
x=linspace(min(SDV(1).Xo),max(SDV(1).Xo),Nnodes);
y=linspace(min(SDV(1).Yo),max(SDV(1).Yo),Nnodes);
[X,Y]=meshgrid(x,y);
rho1=zeros(size(X));
rho2=zeros(size(X));
rho3=zeros(size(X));
rho=zeros(size(X));

for ny=1:length(y)
    for nx=1:length(x)
        rho1(ny,nx)=(feval(SDV(4).fitobject,[X(ny,nx),Y(ny,nx)]));
        rho2(ny,nx)=(feval(SDV(5).fitobject,[X(ny,nx),Y(ny,nx)]));
        rho3(ny,nx)=(feval(SDV(6).fitobject,[X(ny,nx),Y(ny,nx)]));
        
        rho(ny,nx)=(rho1(ny,nx)+rho2(ny,nx)+rho3(ny,nx))./3;
        
    end
    disp(ny);
end

figure(1);
clf;
hold on;
surf(X,Y,rho,'EdgeColor','none');
% csvwrite(['pname','Rho.csv'],rho);

csvwrite([pname '\Rho.csv'],rho);
% csvwrite('Rho.csv',rho);
save RhoIgnore;