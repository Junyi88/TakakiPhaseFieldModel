clear;

Prefix0=['../exec/'];

filename=[Prefix0 'Phi_0.csv'];
Phi=csvread(filename);

filename=[Prefix0 'Theta_0.csv'];
Theta=csvread(filename);

filename=[Prefix0 'Rho_0.csv'];
Rho=csvread(filename);

filename=[Prefix0 'EStored_0.csv'];
EStored=csvread(filename);

figure(1);
clf;
hold on;
surf(Theta,'EdgeColor','none');

figure(2);
clf;
hold on;
surf(Rho,'EdgeColor','none');

figure(3);
clf;
hold on;
surf(EStored,'EdgeColor','none');

%%
filename=[Prefix0 'PhiDx_0.csv'];
PhiDx=csvread(filename);
filename=[Prefix0 'PhiDy_0.csv'];
PhiDy=csvread(filename);
filename=[Prefix0 'PhiDxx_0.csv'];
PhiDxx=csvread(filename);
filename=[Prefix0 'PhiDyy_0.csv'];
PhiDyy=csvread(filename);
filename=[Prefix0 'PhiDxy_0.csv'];
PhiDxy=csvread(filename);
filename=[Prefix0 'PhiD2_0.csv'];
PhiD2=csvread(filename);


filename=[Prefix0 'ThetaDx_0.csv'];
ThetaDx=csvread(filename);
filename=[Prefix0 'ThetaDy_0.csv'];
ThetaDy=csvread(filename);
filename=[Prefix0 'ThetaDxx_0.csv'];
ThetaDxx=csvread(filename);
filename=[Prefix0 'ThetaDyy_0.csv'];
ThetaDyy=csvread(filename);
filename=[Prefix0 'ThetaDxy_0.csv'];
ThetaDxy=csvread(filename);
filename=[Prefix0 'ThetaD2_0.csv'];
ThetaD2=csvread(filename);

%%
filename=[Prefix0 'Wall_dFdPhase_0.csv'];
WdFdPhase=csvread(filename);

filename=[Prefix0 'Bulk_dFdPhase_0.csv'];
BdFdPhase=csvread(filename);

filename=[Prefix0 'Ori_dFdPhase_0.csv'];
OdFdPhase=csvread(filename);

filename=[Prefix0 'Grad_dFdPhase_0.csv'];
GdFdPhase=csvread(filename);

filename=[Prefix0 'dPhidt_0.csv'];
dPhidt=csvread(filename);

filename=[Prefix0 'dThetadt_0.csv'];
dThetadt=csvread(filename);
%%
ypos=100;
figure(4);
clf;
hold on;
plot(Theta(ypos,:),'r-');
plot(ThetaDx(ypos,:),'b-');
plot(ThetaDxx(ypos,:),'k-');
% surf(PhiDx,'EdgeColor','none');



figure(5);
clf;
hold on;
% plot(Theta(ypos,:),'r-');
plot(BdFdPhase(ypos,:),'rx-');
plot(OdFdPhase(ypos,:),'bo-');
plot(GdFdPhase(ypos,:),'k>-');
plot(WdFdPhase(ypos,:),'gx-');

% plot(100.*dThetadt(ypos,:),'k-');