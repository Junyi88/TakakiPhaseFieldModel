function [Estored]=CalculateEstored(Rho,mu,bvec)

    Estored=Rho.*mu.*bvec.*bvec.*0.5;

end