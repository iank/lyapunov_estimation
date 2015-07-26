function [ dxyz ] = lorenz( t, xyz )
    sigma = 10;
    beta = 8/3;
    rho = 28;
    dxyz = zeros(3,1);
    dxyz(1) = sigma*(xyz(2) - xyz(1));
    dxyz(2) = xyz(1)*(rho - xyz(3)) - xyz(2);
    dxyz(3) = xyz(1)*xyz(2) - beta*xyz(3);
end

