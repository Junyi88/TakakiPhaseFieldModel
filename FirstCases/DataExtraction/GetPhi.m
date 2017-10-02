clear;
load names;

ro=csvread('Rho.csv');
Phi=zeros(size(ro));

csvwrite([pname '\Phi.csv'],Phi);
