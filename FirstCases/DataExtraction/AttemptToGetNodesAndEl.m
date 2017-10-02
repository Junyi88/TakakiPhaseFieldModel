clear;

load names;
% NodesRaw.NodeFile='../G49/nodes.inc';
NodesRaw.N=csvread(NodesRaw.NodeFile);


%% Determine How many Z layers

[NodesRaw.Zs]=unique(NodesRaw.N(:,4));
[NNodes,Ndim]=size(NodesRaw.Zs);
Nz=length(NodesRaw.Zs);
Layers=Nz+1;

LayerCount=zeros(Layers,1);

for n1=1:Layers
    Nodes(n1).N=zeros(NNodes,4);
end

for nNode=1:length(NodesRaw.N)
    hit=0;
    for n1=1:Layers-1
        if NodesRaw.N(nNode,4)==NodesRaw.Zs(n1)
            hit=1;
            LayerCount(n1)=LayerCount(n1)+1;
            Nodes(n1).N(LayerCount(n1),:)=NodesRaw.N(nNode,:);
        end
    end
    if hit==0
        LayerCount(Layers)=LayerCount(Layers)+1;
        Nodes(Layers).N(LayerCount(Layers),:)=NodesRaw.N(nNode,:);
    end
end

for n1=1:Layers
    Nodes(n1).N=Nodes(n1).N(1:LayerCount(n1),:);
end

figure(1);
clf;
plot3(NodesRaw.N(:,2),NodesRaw.N(:,3),NodesRaw.N(:,4),'rx')


figure(2);
clf;
hold on;
plot(Nodes(1).N(:,2),Nodes(1).N(:,3),'rx');
plot(Nodes(2).N(:,2),Nodes(2).N(:,3),'bs');

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

%%

figure(200);
clf;
hold on;
nEl=5;
x=zeros(1,20);
y=zeros(1,20);
z=zeros(1,20);

for n1=1:20
    x(n1)=NodesRaw.N(Els(nEl,n1+1),2);
    y(n1)=NodesRaw.N(Els(nEl,n1+1),3);
    z(n1)=NodesRaw.N(Els(nEl,n1+1),4);
end

plot3(x,y,z,'rx-');

plot3(x(1),y(1),z(1),'go-');
plot3(x(end),y(end),z(end),'bo-');

NodePos=NodesRaw.N;
NodeLayers=NodesRaw.Zs;

save NodeDataIgnore NodePos Nodes NodeLayers;