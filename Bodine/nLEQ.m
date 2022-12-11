function dxdt = nLEQ(t,x)
    mu = 3.986004418e5;
    dxdt = [x(2);
    -(mu*x(1))/(x(1)^2+x(3)^2)^(3/2);
    x(4);
    -(mu*x(3))/(x(1)^2+x(3)^2)^(3/2)];
end