% Inputs: x: [1x4][km]
% Note: Only handles for step in time
function [A] = findANominal(x,mu)
if nargin == 1             % Default to the Earth
    mu = 398600;            % [km^3/s^2] Earth Gravitional Parameter
end

A = [0,1,0,0;
    -mu*(x(3)^2-2*x(1)^2)/((x(1)^2+x(3)^2)^(5/2)),0,3*mu*x(1)*x(3)/((x(1)^2+x(3)^2)^(5/2)),0;
    0,0,0,1;
    3*mu*x(3)*x(1)/((x(1)^2+x(3)^2)^(5/2)),0,-mu*(x(1)^2-2*x(3)^2)/((x(1)^2+x(3)^2)^(5/2)),0];
end
