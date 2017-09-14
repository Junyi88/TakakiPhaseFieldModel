function [SeedSize,SeedDist,Seeded,NoTouch,SeedNum]=...
    AddSeedPoint(x,y,SeedSize,SeedDist,Seeded,NoTouch,SeedNum,X,Y)
    
    [NY,NX]=size(Seeded);
    SeedMask=zeros(NY,NX);
    NoTouchMask=zeros(NY,NX);
    NoTouchFlag=0;
    
    for ny=1:NY
        for nx=1:NX
            r=sqrt(((X(ny,nx)-x).^2)+((Y(ny,nx)-y).^2));
            
            if (r<=SeedDist)
               NoTouchFlag=1;
               break;
            elseif (r<=SeedSize)
                SeedMask(ny,nx)=1;
            end
        end
        
        if (NoTouchFlag==1)
            break;
        end
    end
    
    if (NoTouchFlag==1)
       break;
    else
        SeedNum=SeedNum+1;
    end
    
end