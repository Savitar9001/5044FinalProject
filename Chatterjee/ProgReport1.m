% @Author: Sayan Chatterjee
% @Date: 12/1/2022
% @Purpose: The following is code for ASEN5044 Final Project for an
% Orbiting Spacecraft. Divided into System Dynamics and then Kalman Filter
% implementations.
clc; clear; close all;
addpath('Functions');

%% Let's Define Some Constants
mu = 398600;                        % [km^3/s^2] Earth Gravitational Parameters
r0 = 6678;                          % [km] Nominal Orbit of SC
xDot0 = 0;                          % [km/s]
y0 = 0;
yDot0 = r0*sqrt(mu/(r0^3));
rE = 6378;                          % [km] Earth Radius Equatorial
omega = sqrt(mu/r0^3);              % Constant to Simplify Calculation
spinRateEarth = (2*pi)/86400;       % [rad/s] Spinrate of the Earth
detlaT = 10;
numStat = 12;                       % Number of Stations
trans = [0 1 0;-1 0 0; 0 0 1]; % rotation matrix for converting to tangent

% Initial Calculations
theta0 = NaN(1,length(numStat));
for ii = 1:numStat
    theta0(ii) = (ii-1)*(pi/6);
end


%% Linear Simulation
% Forming Time Vector
tVec = linspace(0,14000,1401);



% Find the Nominal Trajectory
xNominal = NaN(length(tVec),4);

for ii = 1:length(tVec)
    xNominal(ii,:) = [r0*cos(omega*tVec(ii)), -omega*r0*sin(omega*tVec(ii)), r0*sin(omega*tVec(ii)), omega*r0*cos(omega*tVec(ii))];
end

% perturbation form solution sketch
pertub = [0.0, 0.075, 0.0, -0.021];
xLin = xNominal(1,:) + pertub;

% Find the Linear Pertub Solution
for kk = 1:(length(tVec)-1)
    A = findANominal(xNominal(kk,:));
    F = eye(4) + detlaT*A;
    xLin(kk+1,:) = xNominal(kk+1,:) + (F*pertub(kk,:)')';
    pertub(kk+1,:) = F*pertub(kk,:)';
end


% Nonlinear Simulation
state0 = xLin(1,:);

% Don't Forget to Up the Tolerance
tol = 1e-12;
options = odeset('RelTol',tol,'AbsTol',[tol tol tol tol]);

[nonLinearT,nonLinearState] = ode45(@propDyDt,tVec,state0,options);





% Iterate Through Time
H = cell(length(tVec),1);                   % Let's Try to Preallocate to be a Good Programmer
xStation = NaN(length(tVec), 2*numStat);
yStation = NaN(length(tVec), 2*numStat);
yNominal = NaN(length(tVec), 3*numStat);
dy = NaN(length(tVec), 3*numStat);
yLinear = NaN(length(tVec), 3*numStat);
yVis = NaN(length(tVec), 3*numStat);
StationNum = NaN(numStat, length(tVec));
yTrue = NaN(length(tVec), 3*numStat);
tTrueVis = NaN(length(tVec), 3*numStat);
trueStation = NaN(numStat, length(tVec));

for ii = 1:length(tVec)
    % Let's Try to Preallocate to be a Good Programmer
    tempX = NaN(1,2*numStat);
    tempY = NaN(1,2*numStat);
    temp_yOut = cell(numStat,1);
    
    
    % generate station measurement data and calculate H matrix
    for kk= 1:numStat
        tempX(2*kk-1:2*kk) = [rE*cos(spinRateEarth*ii*detlaT + theta0(kk));-spinRateEarth*rE*sin(spinRateEarth*ii*detlaT+theta0(kk))];
        tempY(2*kk-1:2*kk) = [rE*sin(spinRateEarth*ii*detlaT + theta0(kk));spinRateEarth*rE*cos(spinRateEarth*ii*detlaT+theta0(kk))];
        
        H{ii} = [H{ii};findHNominal(xNominal(ii,:),tempX(2*kk-1:2*kk),tempY(2*kk-1:2*kk))];

        temp_yOut{kk,1} = [findYMeas(xNominal(ii,:),tempX(2*kk-1:2*kk),tempY(2*kk-1:2*kk))];
    end
    temp_yOut = [temp_yOut{:}]; % Converting Cell to Matrix
   
    
    % store station data
    xStation(ii,:) = tempX;
    yStation(ii,:) = tempY;
    yNominal(ii,:) = temp_yOut(:);
    dy(ii,:) = H{ii}*pertub(ii,:)';

    %check for 'visibility' of satellite  per station
    tempStation = zeros(12,1);
    yLinear(ii,:) = yNominal(ii,:) + dy(ii,:);

    % Quick Preallocation
    tempTheta = NaN(1,numStat);
    tempAngle = NaN(1,numStat);

    for jj =1:numStat
        tempTheta(jj) = atan2(yStation(ii,2*jj-1),xStation(ii,2*jj-1));
        pos = [xStation(ii,2*jj-1),yStation(ii,2*jj-1),0];
        pos2 = [xLin(ii,1) - xStation(ii,2*jj-1),xLin(ii,3)-yStation(ii,2*jj-1),0];
        posh = trans*pos';
        tempAngle(jj) = acos(dot(posh,pos2)/(norm(posh)*norm(pos2)));
        c = cross(posh,pos2);
        yVis(ii,3*jj-2:3*jj) = yLinear(ii,3*jj-2:3*jj);

        if c(3) >= 0
            tempStation(jj) = jj;
        end

        if c(3) < 0
            yVis(ii,3*jj-2:3*jj) = nan;
            tempStation(jj) = nan;
        end

        if yVis(ii,3*jj) < pi
            yVis(ii,3*jj) = yVis(ii,3*jj) + 2*pi;
        end
        if yVis(ii,3*jj) > pi
            yVis(ii,3*jj) = yVis(ii,3*jj) - 2*pi;
        end
    end

    StationNum(:,ii) = tempStation;

end

% Perform Calculations for True Measurements
for ii=1:length(tVec)
    truetemp_y = cell(numStat,1);
    for kk =1:numStat
        truetemp_y{kk} = [findYMeas(nonLinearState(ii,:),xStation(ii,2*kk-1:2*kk),yStation(ii,2*kk-1:2*kk))];
    end
    truetemp_y = [truetemp_y{:}]; % Converting Cell to Matrix

    tempStation = zeros(12,1);
    yTrue(ii,:) = truetemp_y(:);

    for jj =1:12
        tempTheta(jj) = atan2(yStation(ii,2*jj-1),xStation(ii,2*jj-1));
        pos = [xStation(ii,2*jj-1),yStation(ii,2*jj-1),0];
        pos2 = [xLin(ii,1) - xStation(ii,2*jj-1),xLin(ii,3)-yStation(ii,2*jj-1),0];
        posh = trans*pos';
        tempAngle(jj) = acos(dot(posh,pos2)/(norm(posh)*norm(pos2)));
        c = cross(posh,pos2);
        tTrueVis(ii,3*jj-2:3*jj) = yTrue(ii,3*jj-2:3*jj);
        if c(3) >=0
            tempStation(jj) = jj;
        end
        if c(3) < 0
            tempAngle(jj) = -tempAngle(jj);
            tTrueVis(ii,3*jj-2:3*jj) = nan;
            tempStation(jj) = nan;
        end
        if tTrueVis(ii,3*jj) < pi
            tTrueVis(ii,3*jj) = tTrueVis(ii,3*jj) + 2*pi;
        end
        if tTrueVis(ii,3*jj) > pi
            tTrueVis(ii,3*jj) = tTrueVis(ii,3*jj) - 2*pi;
        end
    end

    trueStation(:,ii) = tempStation;

end


% Nonlinear ODE state plot
figure
hold on
grid on
subplot(4,1,1);
plot(tVec,nonLinearState(:,1)); xlabel('Time (seconds)'); ylabel('X [km]');
title('States vs Time, Full Nonlinear Dynamics Simulation');
subplot(4,1,2);
plot(tVec,nonLinearState(:,2)); xlabel('Time (seconds)'); ylabel('$\dot{X} [km/s]$','Interpreter','latex');
subplot(4,1,3);
plot(tVec,nonLinearState(:,3)); xlabel('Time (seconds)'); ylabel('Y [km]');
subplot(4,1,4);
plot(tVec,nonLinearState(:,4)); xlabel('Time (seconds)'); ylabel('$\dot{Y} [km/s]$','Interpreter','latex');
%



figure
hold on
grid on
subplot(4,1,1);
plot(tVec,xLin(:,1)); xlabel('Time (seconds)'); ylabel('X [km]'); 
title('States vs Time, Linearized Approximate Dynamics Simulation');
subplot(4,1,2);
plot(tVec,xLin(:,2)); xlabel('Time (seconds)'); ylabel('$\dot{X} [km/s]$','Interpreter','latex');
subplot(4,1,3);
plot(tVec,xLin(:,3)); xlabel('Time (seconds)'); ylabel('Y [km]');
subplot(4,1,4);
plot(tVec,xLin(:,4)); xlabel('Time (seconds)'); ylabel('$\dot{Y} [km/s]$','Interpreter','latex');



figure
hold on
grid on
subplot(4,1,1);
plot(tVec,xLin(:,1)-xNominal(:,1)); xlabel('Time (seconds)'); ylabel('\delta X [km]');
title('Linearized Approximate Perturbations vs Time');
subplot(4,1,2);
plot(tVec,xLin(:,2)-xNominal(:,2)); xlabel('Time (seconds)'); ylabel('$\dot{X} [km/s]$','Interpreter','latex');
subplot(4,1,3);
plot(tVec,xLin(:,3)-xNominal(:,3)); xlabel('Time (seconds)'); ylabel(' \delta Y [km]'); 
subplot(4,1,4);
plot(tVec,xLin(:,4)-xNominal(:,4)); xlabel('Time (seconds)'); ylabel('$\dot{Y} [km/s]$','Interpreter','latex');






figure;
subplot(4,1,1);

for ii = 1:12
    plot(tVec,tTrueVis(:,3*ii-2), 'x', 'markersize', 4); hold on; grid on; ylabel('\rho^i [km]'); 
end
title('Full Nonlinear Model Data Simulation');
subplot(4,1,2);
for ii = 1:12
    plot(tVec,tTrueVis(:,3*ii-1), 'o', 'markersize', 4); hold on; grid on; ylabel('$\dot{\rho^i} [km/s]$','Interpreter','latex');
end

subplot(4,1,3);
for ii = 1:12
    plot(tVec,tTrueVis(:,3*ii), 'o', 'markersize', 4); hold on; grid on; ylabel(' \phi ^i [rads]');
end
subplot(4,1,4);
for ii = 1:12
    plot(tVec,trueStation(ii,:), '.'); hold on; ylabel('Visible Station ID'); 
end





figure;
subplot(4,1,1);
for ii = 1:12
    plot(tVec,yVis(:,3*ii-2), 'x', 'markersize', 4); hold on; grid on; ylabel('\rho^i [km]'); 
end
title('Approximate Linearized Model Data Simulation');
subplot(4,1,2);
for ii = 1:12
    plot(tVec,yVis(:,3*ii-1), 'o', 'markersize', 4); hold on; grid on; ylabel('$\dot{\rho^i} [km/s]$','Interpreter','latex');
end

subplot(4,1,3);
for ii = 1:12
    plot(tVec,yVis(:,3*ii), 'o', 'markersize', 4); hold on; grid on; ylabel(' \phi ^i [rads]');
end

subplot(4,1,4);
for ii = 1:12
    plot(tVec,StationNum(ii,:),'.'), hold on; grid on; ylabel('Visible Station ID'); 
end



