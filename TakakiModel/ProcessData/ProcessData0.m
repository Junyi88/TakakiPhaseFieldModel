clear;

Prefix0=['../exec/'];

filename=[Prefix0 'Phi_0.csv'];
Phi=csvread(filename);

filename=[Prefix0 'Theta_0.csv'];
Theta=csvread(filename);

filename=[Prefix0 'Rho_0.csv'];
Rho=csvread(filename);

figure(1);
clf;
hold on;
surf(Rho,'EdgeColor','none');