function [Seeded,SeedList,SeedNum]=...
    GenerateSeedList(EVal,X,Y,EMin)

    [NY,NX]=size(EVal);
    SeedList=zeros(NY.*NX,3);
    
    xList=X(:);
    yList=Y(:);
    EList=EVal(:);
    
    Seeded=zeros(NY,NX);
    SeedNum=0;
    LList=length(EList);
    
    for n1=1:LList
        if (EList(n1)>=EMin)
            SeedNum=SeedNum+1;
            SeedList(SeedNum,1)=xList(n1);
            SeedList(SeedNum,2)=yList(n1);
            SeedList(SeedNum,3)=EList(n1);
        end
    end
    
    %%
    SeedList=SeedList(1:SeedNum,:);
    [~,I] = sort(SeedList(:,3),'descend');
    
    SeedList=SeedList(I,:);
    SeedNum=0;
end