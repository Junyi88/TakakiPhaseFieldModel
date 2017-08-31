clear;

ElRaw.ElFile='../G49/elements.inc';
ElRaw.ElRaw=csvread(ElRaw.ElFile);
[Ny,Nx]=size(ElRaw.ElRaw);

NyReal=Ny/2;
NxReal=Nx-2;

Els=zeros(NyReal,NxReal+1);

count=1;
for n1=1:NyReal
    Els(n1,1:16)=ElRaw.ElRaw(count,1:16);
    count=count+1;
    
    Els(n1,17:21)=ElRaw.ElRaw(count,1:5);
    count=count+1;
end
