clear;

Filename='../G49/rotationMatrix.inc';

fid = fopen(Filename);
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

step=0;
EulerAngles=zeros(length(g),3);
count=0;

TF = strcmp(g(2).s(2:9),'MATERIAL');

for n1=1:length(g)
    
    if length(g(n1).s)>=9
    if (step==0)&&(strcmp(g(n1).s(2:9),'MATERIAL'))
        count=count+1;
        step=1;
    elseif (step==1)
        step=0;
        EulerAngles(count,:)=str2num(g(n1).s(3:end));
    end
    end
    
end

EulerAngles=EulerAngles(1:count,:);

save EulerAnglesIgnore EulerAngles;