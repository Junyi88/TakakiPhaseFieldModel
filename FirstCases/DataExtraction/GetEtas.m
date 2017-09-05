clear;
ro=csvread('Rho.csv');
Eta=zeros(size(ro));

csvwrite('Eta.csv',Eta);