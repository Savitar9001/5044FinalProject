function [y] = measureY(x,XStation,YStation)
% calculates y measurements from given equations
% x is state, XStation, YStation is station location
y = [sqrt((x(1)-XStation(1))^2 + (x(3)-YStation(1))^2);
    ((x(1) - XStation(1))*(x(2) - XStation(2))+(x(3) - YStation(1))*(x(4)-YStation(2)))/sqrt((x(1)-XStation(1))^2 + (x(3)-YStation(1))^2);
    atan2((x(3)-YStation(1)),(x(1)-XStation(1)))];
end
