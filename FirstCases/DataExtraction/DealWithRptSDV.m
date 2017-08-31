clear;

SDVFilename='../G49/SDVNodal.rpt';

fid = fopen(SDVFilename);
count=1;
tline = fgetl(fid);
g(count).s=tline;
while ischar(tline)
%     disp(tline)
    tline = fgetl(fid);
    count=count+1;
    g(count).s=tline;
end

fclose(fid);

%%

NLines=length(g);
SDVS=zeros(NLines,4);
count=0;
for n1=1:NLines-1
    X=str2num(g(n1).s);
    if length(X)==4
        count=count+1;
        SDVS(count,:)=X;
    end
    
end

SDVS2=SDVS(1:count,:);
G=unique(SDVS2(:,1));

NodeCounts=zeros(length(G),1);
Fx=zeros(length(G),4);
for n1=1:length(SDVS2)
   NodeCounts(SDVS2(n1,1))= NodeCounts(SDVS2(n1,1))+1;
   Fx(SDVS2(n1,1),:)=SDVS2(n1,:);
end

figure(1);
clf;
plot(NodeCounts,'rx-');

NN=Fx(8749,:,:);

SDVs=Fx;
save SDVDataIgnore SDVs;