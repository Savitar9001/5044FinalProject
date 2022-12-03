function [H] = findHNominal(x,XS,YS)
% calculates H matrix  jacobian
% x is state, XSS, YSS is station location
y1 = sqrt((x(1)-XS(1))^2 + (x(3)-YS(1))^2);
H = [(x(1)-XS(1))/y1,0,(x(3)-YS(1))/y1,0;
    (x(2)-XS(2))/y1-((x(1) - XS(1))*(x(2)-XS(2))+(x(3)-YS(1))*(x(4)-YS(2)))*(x(1)- XS(1))/(y1^3),(x(1)-XS(1))/y1,...
    (x(4)-YS(2))/y1-((x(3) - YS(1))*(x(4)-YS(2))+(x(1) - XS(1))*(x(2)-XS(2)))*(x(3)- YS(1))/(y1^3),(x(3)-YS(1))/y1;
    -(x(3)-YS(1))/((x(3) - YS(1))^2 + (x(1)-XS(1))^2),0,(x(1)-XS(1))/((x(3)-YS(1))^2+(x(1)- XS(1))^2),0];
end

