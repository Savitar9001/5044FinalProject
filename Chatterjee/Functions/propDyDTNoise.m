function xdot_nonlin = propDyDTNoise(t,stateNonLinear)

global wk                       % I hate Using This, but Here is the Noise! 

mu = 3.98600e5;                 % [km^3/s^2] Earth Gravitional Parameter
x1 = stateNonLinear(1);
x2 = stateNonLinear(2);
x3 = stateNonLinear(3);
x4 = stateNonLinear(4);

r = norm([x1 x3]);

xdot_nonlin = [x2; -mu*x1/r^3 + wk(1); x4; -mu*x3/r^3 + wk(2)];
end