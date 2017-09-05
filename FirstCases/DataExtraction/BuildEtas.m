clear;

Eta0=csvread('Eta.csv');

Eta0p1=0.1.*ones(size(Eta0));
Eta0p2=0.1.*ones(size(Eta0));
Eta1=1.0.*ones(size(Eta0));

csvwrite('Eta0p1.csv',Eta0p1);
csvwrite('Eta0p2.csv',Eta0p2);
csvwrite('Eta1.csv',Eta1);