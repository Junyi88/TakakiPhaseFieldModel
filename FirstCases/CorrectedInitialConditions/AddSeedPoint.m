function [Seeded,SeedList,SeedNum]=...
    AddSeedPoint(SeedSize,SeedDist,Seeded,SeedList,SeedNum,X,Y)
    
    [NY,NX]=size(Seeded);
    SeedNum=SeedNum+1;
    xSeed=SeedList(1,1);
    ySeed=SeedList(1,2);
    
    %%
    for ny=1:NY
       for nx=1:NX
           x=X(ny,nx);
           y=Y(ny,nx);
           
           r=sqrt((x-xSeed).^2+(y-ySeed).^2);
           
           if (r<=SeedSize)
               Seeded(ny,nx)=SeedNum;
           end
       end
    end
    
    %%
    [L,~]=size(SeedList);
    Keep=true(L,1);
    Keep(1)=0;
    for n1=2:L
         x=SeedList(n1,1);
         y=SeedList(n1,2);
         r=sqrt((x-xSeed).^2+(y-ySeed).^2);
         
         if (r<=SeedDist)
            Keep(n1,1)=false;
         end
    end
    
    %%
    SeedList=SeedList(Keep,:);
end

