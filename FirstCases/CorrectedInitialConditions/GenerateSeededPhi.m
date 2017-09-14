%% Setup Parameters (User Input)
clear;

Theta=csvread('EulerAngles.csv'); % Euler Angle
Rho=csvread('Rho.csv'); % Rho

mu=25.0e3; % Shear Modulus
bvec=0.28; % Burgess Vector

Ecritical=0.7; % Critical Stored Energy to start Seed
DThetaCrit=0.1;  % Minimum gradient of euler angle to seed
EMin=1e-8; % Minimum Stored Energy to be swtiched on for seeding (set to a small value)

SeedSize=3; % radius of seeds
SeedDist=10; % Minimum Distance Between seeds

PhiMin=0.0; % Phi for not seeds
PhiMax=1.0; % Phi for Seeds

h=1.0; % Grid Size 
%% Initialise Matrix
[NY,NX]=size(Theta); 
x=1:h:NX;
y=1:h:NY;

[X,Y]=meshgrid(x,y);

%% Calculate Stored Energies and Magnitude of misorientation 
[Estored]=CalculateEstored(Rho,mu,bvec);
sTheta=fit([X(:) Y(:)], Theta(:),'linearinterp');
[ThetaX,ThetaY,ThetaXX,ThetaXY,ThetaYY]=differentiate(sTheta,X,Y);
DTheta=sqrt(ThetaX.^2+ThetaY.^2); % misorientation 

%% Calculate 
ELarge=zeros(NY,NX); % 1 if Estored is larger than Ecritical
ThetaLarge=zeros(NY,NX); % 1 if Misorientation is larger than critical Misorientation

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

EThetaLarge=ELarge.*ThetaLarge; % 1 if both Estored and misorientation is larger than critical values
ESeeds=EThetaLarge.*Estored; % EStored values of seeded points
EVal=ESeeds; % Same as previous

%% Perform Seeding
[Seeded,SeedList,SeedNum]=GenerateSeedList(EVal,X,Y,EMin);
[ListLength,~]=size(SeedList);

while (ListLength>0)

[Seeded,SeedList,SeedNum]=...
    AddSeedPoint(SeedSize,SeedDist,Seeded,SeedList,SeedNum,X,Y);
[ListLength,~]=size(SeedList);
end

% Seeded is the variable which recodes the seed number at each location

%% Calculate Phi
Phi=zeros(NY,NX);
for ny=1:NY
    for nx=1:NX
        if (Seeded(ny,nx)>0.5)
            Phi(ny,nx)=PhiMax;
        else
            Phi(ny,nx)=PhiMin;
        end
    end
end

csvwrite('PhiSeeded.csv',Phi);

%% Plot Curves
figure(1);
clf;
subplot(2,2,1);
surf(X,Y,Estored,'EdgeColor','none');
title(['\phi ']);
view([0 90]);
xlim([0 211]);
ylim([0 211]);
pbaspect([1 1 1]);


subplot(2,2,2);
surf(X,Y,DTheta,'EdgeColor','none');
title(['|\nabla \theta|']);
view([0 90]);
xlim([0 211]);
ylim([0 211]);
pbaspect([1 1 1]);

subplot(2,2,3);
surf(X,Y,EVal,'EdgeColor','none');
title(['EStored of Possible Seeding Locations']);
view([0 90]);
xlim([0 211]);
ylim([0 211]);
pbaspect([1 1 1]);

subplot(2,2,4);
surf(X,Y,Seeded,'EdgeColor','none');
title(['Seed Number']);
view([0 90]);
xlim([0 211]);
ylim([0 211]);
pbaspect([1 1 1]);

figure(2);
clf;
surf(X,Y,Phi,'EdgeColor','none');
title(['\phi']);
view([0 90]);
xlim([0 211]);
ylim([0 211]);
pbaspect([1 1 1]);