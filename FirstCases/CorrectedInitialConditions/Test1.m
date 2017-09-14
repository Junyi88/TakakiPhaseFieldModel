clear;

Theta=csvread('EulerAngles.csv');
Rho=csvread('Rho.csv');

mu=25.0e3;
bvec=0.28;

[NY,NX]=size(Theta);
x=1:NX;
y=1:NY;

[X,Y]=meshgrid(x,y);
Ecritical=0.7;
DThetaCrit=1.0e-3;
%%
[Estored]=CalculateEstored(Rho,mu,bvec);

figure(1);
clf;
hold on;
surf(Estored,'EdgeColor','none');
pbaspect([1 1 1]);

figure(2);
clf;
hold on;
surf(Theta,'EdgeColor','none');
pbaspect([1 1 1]);

sTheta=fit([X(:) Y(:)], Theta(:),'linearinterp');

[ThetaX,ThetaY,ThetaXX,ThetaXY,ThetaYY]=differentiate(sTheta,X,Y);

figure(3);
clf;
hold on;
surf(ThetaX,'EdgeColor','none');
pbaspect([1 1 1]);



%% Calculate differentiation of YY

DTheta=sqrt(ThetaX.^2+ThetaY.^2);

figure(3);
clf;
hold on;
surf(DTheta,'EdgeColor','none');
pbaspect([1 1 1]);

figure(8);
clf;
hold on;
hist(Estored(:));
pbaspect([1 1 1]);

ELarge=zeros(NY,NX);
ThetaLarge=zeros(NY,NX);


for ny=1:NY
    for nx=1:NX
        if (Estored(ny,nx)>=Ecritical)
            ELarge(ny,nx)=1;
        end
        if (DTheta(ny,nx)>=DThetaCrit)
            ThetaLarge(ny,nx)=1;
        end
    end
end
EThetaLarge=ELarge.*ThetaLarge;



ESeeds=EThetaLarge.*Estored;

figure(4);
clf;
hold on;
surf(ESeeds,'EdgeColor','none');
pbaspect([1 1 1]);