function xdot_nonlin = ode_nonlin(t,state_nonlin)
    mu = 3.98600e5;
    x1 = state_nonlin(1);
    x2 = state_nonlin(2);
    x3 = state_nonlin(3);
    x4 = state_nonlin(4);
    
    r = norm([x1 x3]);
    
    xdot_nonlin = [x2; -mu*x1/r^3; x4; -mu*x3/r^3];
end