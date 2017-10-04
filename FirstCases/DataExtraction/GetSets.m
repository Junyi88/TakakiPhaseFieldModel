%% 2/9 to run: Get sets.inc

clear;
load names;
Filename=GetSetsFilename;
% Filename='../G49/sets.inc';

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

count=0;
step=0;
for n1=1:length(g)
    if length(g(n1).s)>=13
        if (step==0)&&(strcmp([g(n1).s(1:6) g(n1).s(8:13)],'*Nset,nset=g'))
            count=count+1;
            step=1;
        elseif (step==1)
            
            SetNodes(count).N=str2num(g(n1).s(1:end));
            step=2;
        elseif (step==2)
            SetNodes(count).N=[SetNodes(count).N str2num(g(n1).s(1:end))];
            if (strcmp(g(n1).s(1),'*'))
                step=0;
                SetNodes(count).N=SetNodes(count).N(SetNodes(count).N~=0);
            end
        end
    end
end

save SetsIgnore SetNodes;