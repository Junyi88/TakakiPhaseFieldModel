%% 8/9 to run: Write phi in .csv
clear;
load names;
load cut;

ro=csvread('Rho.csv');
Phi=zeros(size(ro));

% Phi=Phi(Ysmall+1:Ylarge-1,Xsmall+1:Xlarge-1);

csvwrite([pname '\Phi.csv'],Phi);
