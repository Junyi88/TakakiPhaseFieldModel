%% 7/9 to run: Write rho in .csv
clear;

load SetsIgnore;
load EulerAnglesIgnore;
load SDVSortedIgnore;
load NodeDataIgnore;

%%
load names;

Nnodes=numberOfNodes;


x=min(SDV(5).X):dx:max(SDV(5).X);
y=min(SDV(5).Y):dx:max(SDV(5).Y);

[X,Y]=meshgrid(x,y);
rho1=zeros(size(X));
rho2=zeros(size(X));
rho3=zeros(size(X));
rho=zeros(size(X));

% %Cut NoN number
Xsmall=1+0; %Cut 0 column
Xlarge=length(x)-0;
Ysmall=1+0;
Ylarge=length(y)-0;

for ny=1:length(y)
    for nx=1:length(x)
        rho1(ny,nx)=(feval(SDV(4).fitobject,[X(ny,nx),Y(ny,nx)]));
        rho2(ny,nx)=(feval(SDV(5).fitobject,[X(ny,nx),Y(ny,nx)]));
        rho3(ny,nx)=(feval(SDV(6).fitobject,[X(ny,nx),Y(ny,nx)]));
        
        rho(ny,nx)=(rho1(ny,nx)+rho2(ny,nx)+rho3(ny,nx))./3;
    end
    disp(ny);
end

X=X(Ysmall+1:Ylarge-1,Xsmall+1:Xlarge-1);
Y=Y(Ysmall+1:Ylarge-1,Xsmall+1:Xlarge-1);
rho=rho(Ysmall+1:Ylarge-1,Xsmall+1:Xlarge-1);



figure(71);
clf;
hold on;
surf(X,Y,rho,'EdgeColor','none');

csvwrite([pname '\Rho.csv'],rho);
% csvwrite('Rho.csv',rho);
save RhoIgnore;
save cut Xsmall Xlarge Ysmall Ylarge;