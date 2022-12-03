function DyDt = propDyDt(t,state)
    mu = 398600;            % [km^3/s^2] Earth Gravitional Parameter
    x1 = state(1);
    x2 = state(2);
    x3 = state(3);
    x4 = state(4);
    
    r = norm([x1 x3]);
    
    DyDt = [x2; -mu*x1/r^3; x4; -mu*x3/r^3];
end