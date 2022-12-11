function out = A(x,mu)
    r2 = x(1)^2 + x(3)^2; % This is r^2
    
    f2x1 = -mu*((r2^(3/2))-(3*x(1)^2*sqrt(r2)))/(r2^3);
    f4x3 = -mu*((r2^(3/2))-(3*x(3)^2*sqrt(r2)))/(r2^3);
    f4x1 = (3*mu*x(3)*x(1))/(r2^(5/2));
    f2x3 = (3*mu*x(3)*x(1))/(r2^(5/2));
    
    out = [0   1   0  0
        f2x1 0 f2x3 0;
         0   0   0  1;
        f4x1 0 f4x3 0];
end