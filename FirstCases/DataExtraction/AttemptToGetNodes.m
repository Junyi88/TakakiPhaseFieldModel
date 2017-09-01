clear;

NodesRaw.NodeFile='../G49/nodes.inc';
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