clear;

F=csvread('Test1_F.csv');

figure(1);
clf;
hold on;
surf(F,'EdgeColor','none');

Dx=csvread('Test1_Dx.csv');
Dy=csvread('Test1_Dy.csv');
Dxx=csvread('Test1_Dxx.csv');
Dyy=csvread('Test1_Dyy.csv');
Dxy=csvread('Test1_Dxy.csv');
D2=csvread('Test1_D2.csv');

Dxa=csvread('../data/Cor0_Fx.csv');
Dya=csvread('../data/Cor0_Fy.csv');
Dxxa=csvread('../data/Cor0_Fxx.csv');
Dyya=csvread('../data/Cor0_Fyy.csv');
Dxya=csvread('../data/Cor0_Fxy.csv');
D2a=csvread('../data/Cor0_F2.csv');

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