clear;

F=csvread('Test1B_F.csv');

figure(1);
clf;
hold on;
surf(F,'EdgeColor','none');

Dx=csvread('Test1B_Dx.csv');
Dy=csvread('Test1B_Dy.csv');
Dxx=csvread('Test1B_Dxx.csv');
Dyy=csvread('Test1B_Dyy.csv');
Dxy=csvread('Test1B_Dxy.csv');
D2=csvread('Test1B_D2.csv');

Dxa=csvread('../data/Cor0p005_Fx.csv');
Dya=csvread('../data/Cor0p005_Fy.csv');
Dxxa=csvread('../data/Cor0p005_Fxx.csv');
Dyya=csvread('../data/Cor0p005_Fyy.csv');
Dxya=csvread('../data/Cor0p005_Fxy.csv');
D2a=csvread('../data/Cor0p005_F2.csv');

x=linspace(1,201,201);
y=linspace(1,201,201);
[X,Y]=meshgrid(x,y);

figure(2);
clf;
hold on;
surf(Dx,'EdgeColor','none');
plot3(X(:),Y(:),Dxa(:),'rx');

figure(3);
clf;
hold on;
surf(Dy,'EdgeColor','none');
plot3(X(:),Y(:),Dya(:),'rx');

figure(4);
clf;
hold on;
surf(Dxx,'EdgeColor','none');
plot3(X(:),Y(:),Dxxa(:),'rx');

figure(5);
clf;
hold on;
surf(Dyy,'EdgeColor','none');
plot3(X(:),Y(:),Dyya(:),'rx');

figure(6);
clf;
hold on;
surf(Dxy,'EdgeColor','none');
plot3(X(:),Y(:),Dxya(:),'rx');

figure(7);
clf;
hold on;
surf(D2,'EdgeColor','none');
plot3(X(:),Y(:),D2a(:),'rx');