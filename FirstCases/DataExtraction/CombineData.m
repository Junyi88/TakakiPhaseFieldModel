clear;

load SetsIgnore;
load EulerAnglesIgnore;
load SDVSortedIgnore;
load NodeDataIgnore;

%%
Nnodes=211;
x=linspace(min(SDV(1).Xo),max(SDV(1).Xo),Nnodes);
y=linspace(min(SDV(1).Yo),max(SDV(1).Yo),Nnodes);
[X,Y]=meshgrid(x,y);
S=zeros(size(X));
A=zeros(size(X));
NodeSet1.x=NodePos(:,2);
NodeSet1.y=NodePos(:,3);
NodeSet1.s=zeros(size(NodeSet1.y));

for n1=1:length(SetNodes)
   for n2=1:length(SetNodes(n1).N)
       n=SetNodes(n1).N(n2);
       NodeSet1.s(n)=n1;
   end
end

No5=Nodes(5).N(:,1);
x5=Nodes(5).N(:,2);
y5=(Nodes(5).N(:,3));
s5=zeros(size(No5));
angle5=zeros(size(No5));
for n1=1:length(No5)
   s5(n1)=NodeSet1.s(No5(n1));
   
   temp=pi.*EulerAngles(s5(n1),3)./180;
   angle5(n1)=wrapToPi(temp);
   
end

Surf1=fit([x5 y5],s5,'linearinterp');
Surf2=fit([x5 y5],angle5,'linearinterp');
for ny=1:length(y)
    for nx=1:length(x)
        S(ny,nx)=(feval(Surf1,[X(ny,nx),Y(ny,nx)]));
        A(ny,nx)=(feval(Surf2,[X(ny,nx),Y(ny,nx)]));
    end
    disp(ny);
end

S2=S;

for ny=1:length(y)
    for nx=1:length(x)
        f1=abs(S(ny,nx)-floor(S(ny,nx)));
        f2=abs(S(ny,nx)-ceil(S(ny,nx)));
        if (f1>0.001)&&(f2>0.001)
            S2(ny,nx)=floor(S(ny,nx));
        end
    end
end

figure(1);
clf;
hold on;
surf(X,Y,S2,'EdgeColor','none')

figure(2);
clf;
hold on;
surf(X,Y,A,'EdgeColor','none');

save EulerAnglesSorted1;