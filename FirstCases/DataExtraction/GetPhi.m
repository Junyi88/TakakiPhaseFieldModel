%% 8/9 to run: Write phi in .csv
clear;
load names;

ro=csvread('Rho.csv');
Phi=zeros(size(ro));

csvwrite([pname '\Phi.csv'],Phi);
