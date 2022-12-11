function out = C(x,t,i)
    Re = 6378; % [km]
    we = 2*pi/86400; % [rad/s]
    thi = (i - 1)*pi/6; % Theta 0 for ground station i
    Xis = Re*cos(we*t+thi);
    Yis = Re*sin(we*t+thi);
    
    dXis = -Re*we*sin(we*t+thi);
    dYis = Re*we*cos(we*t+thi);
    
    sq = sqrt((x(1) - Xis)^2 + (x(3)-Yis)^2);
    
    h1x1 = (x(1)-Xis)/sq;
    h1x3 = (x(3)-Yis)/sq;
    h2x1 = (x(2)-dXis)/sq - (x(1)-Xis)*(x(2)-dXis)/sq^3 - (x(1)-Xis)*(x(3)-Yis)*(x(4)-dYis)/sq^3;
    h2x2 = (x(1)-Xis)/sq;
    h2x3 = (x(4)-dYis)/sq - (x(3)-Yis)*(x(4)-dYis)/sq^3 - (x(1)-Xis)*(x(3)-Yis)*(x(2)-dXis)/sq^3;
    h2x4 = (x(3)-Yis)/sq;
    h3x1 = -(x(3)-Yis)/sq^2;
    h3x3 = (x(1)-Xis)/sq^2;
    
    out = [h1x1 0 h1x3 0;
           h2x1 h2x2 h2x3 h2x4;
           h3x1 0 h3x3 0];
end