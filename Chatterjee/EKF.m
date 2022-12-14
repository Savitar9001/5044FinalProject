% @Author: Sayan Chatterjee
% @Date: 12/1/2022
% @Purpose: The following is code for ASEN5044 Final Project for an
% Orbiting Spacecraft. Divided into System Dynamics and then Kalman Filter
% implementations.
% @Note: The Following is the Linear Kalman Filter Implementation
clc; clear; close all;
addpath('Functions');

%% Load Data
load('orbitdeterm_finalproj_KFdata.mat');

global wk; % Yes, Yes. I know. I thought it was the best idea of using global for noise.
rng(100);

%% Initial Parameters
t0              = 0;
state0          = 6678;                             % [km] Nominal Orbit of SC
y0              = 0;
mu              = 398600;                           % [km^3/s^2] Earth Gravitational Parameters
r0              = sqrt(state0^2 + y0^2);
yDot0           = r0*sqrt(mu/(r0^3));
rE              = 6378;                             % [km] Earth Radius Equatorial
omega           = sqrt(mu/r0^3);
spinRateEarth   = (2*pi)/86400;                     % [rad/s] Spinrate of the Earth
deltaT          = 10;
dtVec           = [0,10];
Per             = ceil(2*pi*sqrt(r0^3/mu));
steps           = 1401;
tVec            = linspace(0,14000,steps);          % Time Vector
TA              = linspace(0,2*pi,steps);           % [rad] True Anamoly
stationNum      = 12;
trans           = [0 1 0;-1 0 0; 0 0 1];
simsN           = 50;


%% Start Prelminary Calculations
theta0 = NaN(stationNum,1);     % Preallocate to Speed Up
for ii = 1:stationNum
    theta0(ii) = (ii-1)*(pi/6);                     % Calculate Distribution of the Stations
end

gamma = [0 0; 1 0; 0 0; 0 1];                       % Disturbance model

% EKF matrices
QTune       = Qtrue;
omegaEXT    = deltaT*gamma;

dx          = [0.01, 0.075, 0.01, -0.021];

% Initialize for KF
xNom    = [r0*cos(omega*t0), -omega*r0*sin(omega*t0), r0*sin(omega*t0), omega*r0*cos(omega*t0)];
P0      = diag((dx./2).^2);
% dx_rand = mvnrnd(zeros(1,4),P0);
dx_rand = [0, 0, 0, 0]; % only use this line for part 3
x0      = [r0*cos(omega*t0), -omega*r0*sin(omega*t0), r0*sin(omega*t0), omega*r0*cos(omega*t0)] + dx_rand;

% Configure ODE Really Quick
tol = 1e-12;
options = odeset('RelTol',tol,'AbsTol',[tol tol tol tol]);


%% Iteration of Sims
for n = 1:simsN
    xTrue = NaN(steps,size(x0,2));
    xTrue(1,:) = x0;                          % Initial estimate

    % simulation with process noise
    for ii = 2:length(tVec)
        wk = mvnrnd(zeros(1,2),Qtrue);
        [t_nonlin,x_noise] = ode45(@propDyDTNoise,dtVec,xTrue(ii-1,:),options);
        xTrue(ii,:) = x_noise(end,:);
    end

    xP = zeros(length(tVec),4);
    xHatm = zeros(length(tVec)-1,4);
    xP(1,:) = xNom;
    P_predict = P0;
    stDev = 2*sqrt(diag(P_predict))';


    %% Let's Iterate Through Time
    H                   = cell(length(tVec)-1,1);                   % Let's Try to Preallocate to be a Good Programmer
    yStackNoise         = NaN(length(tVec), 3*stationNum);
    xStation            = NaN(length(tVec), 2*stationNum);
    yStation            = NaN(length(tVec), 2*stationNum);
    StationNumber       = NaN(stationNum,length(tVec));
    yStoreNominal       = cell(length(tVec),1);
    yStoreNoise         = cell(length(tVec),1);
    errorY              = cell(length(tVec)-1,1);
    NEES                = zeros(length(tVec),1);
    NIS                 = zeros(length(tVec),1);

    for ii = 2:length(tVec)
        [t_nonlin,xNext]   = ode45(@ode_nonlin,dtVec,xP(ii-1,:),options);
        xHatm(ii,:)         = xNext(end,:);
        A                   = findANominal(xP(ii-1,:));                   % CT system matrix at each time step
        F                   = eye(4) + deltaT*A;                % DT system matrix at each time step
        phatm               = (F*P_predict*F' + omegaEXT*QTune*omegaEXT');

        tempX           = NaN(1,2*stationNum);
        tempY           = NaN(1,2*stationNum);
        H{ii}           = [];
        temp_yOut       = cell(stationNum,1);
        temp_yNoise     = cell(stationNum,1);

        for kk = 1:stationNum
            tempX(2*kk-1:2*kk)  = [rE*cos(spinRateEarth*ii*deltaT + theta0(kk));-spinRateEarth*rE*sin(spinRateEarth*ii*deltaT+theta0(kk))];
            tempY(2*kk-1:2*kk)  = [rE*sin(spinRateEarth*ii*deltaT + theta0(kk));spinRateEarth*rE*cos(spinRateEarth*ii*deltaT+theta0(kk))];

            temp_yOut{kk,1}     = [measureY(xHatm(ii,:),tempX(2*kk-1:2*kk),tempY(2*kk-1:2*kk))];
            temp_yNoise{kk,1}   = [measureY(xTrue(ii,:),tempX(2*kk-1:2*kk),tempY(2*kk-1:2*kk))+mvnrnd(zeros(1,3),Rtrue)'];
        end
        temp_yOut   = [temp_yOut{:}]; % Converting Cell to Matrix
        temp_yNoise = [temp_yNoise{:}];

        xStation(ii,:)      = tempX;
        yStation(ii,:)      = tempY;
        yNominal            = temp_yOut(:);
        y_noise             = temp_yNoise(:);
        yStackNoise(ii,:)   = yNominal;

        tempStation = zeros(12,1);

        for j =1:stationNum
            pos = [xStation(ii,2*j-1),yStation(ii,2*j-1),0];                                    % station position
            pos2 = [xTrue(ii,1) - xStation(ii,2*j-1),xTrue(ii,3)-yStation(ii,2*j-1),0];         % satellite position
            posh = trans*pos';                                                                  % position vector rotated to make it tangent
            c = cross(posh,pos2);

            if c(3) >= 0
                tempStation(j) = j;
                H{ii} = [H{ii};findHNominal(xHatm(ii,:),tempX(2*j-1:2*j),tempY(2*j-1:2*j))];
            end

            if c(3) < 0
                yNominal(3*j-2:3*j) = 0;
                y_noise(3*j-2:3*j) = 0;
            end

            yNominal(3*j) = wrapToPi(yNominal(3*j));
            y_noise(3*j) = wrapToPi(y_noise(3*j));

        end

        StationNumber(:,ii) = tempStation;
        yStoreNominal{ii} = yNominal(find(yNominal~=0));
        yStoreNoise{ii} = y_noise(find(y_noise~=0));

        %% Implementation of the EKF

        if length(yStoreNoise{ii}) > 1 % make sure we have data
            %        if length(ydata{ii}) > 1 % make sure we have data

            errorY{ii} = (yStoreNoise{ii} - yStoreNominal{ii});
            if size(yStoreNominal{ii},1) > 3 && size(ydata{ii}(1:3,:),2) > 1
                currYData = ydata{ii}(1:3,:); currYData = currYData(:);
            elseif size(yStoreNominal{ii},1) <= 3
                currYData = ydata{ii}(1:3,end); currYData = currYData(:);
            elseif size(yStoreNominal{ii},1) > 3 || size(ydata{ii}(1:3,:),2) ~= 2
                currYData = ydata{ii}(1:3,1); currYData = currYData(:);
                % yStoreNominal{ii} = yStoreNominal{ii}(1:3);
                currYData = [currYData;0;0;0];
            end
            errorY{ii}  = currYData - yStoreNominal{ii};

            for z = 1:length(errorY{ii})/3 % wrap per measurement
                temp_err_y = errorY{ii};
                temp_err_y(3*z) = wrapToPi(temp_err_y(3*z));
                errorY{ii} = temp_err_y;
            end

            % Compute EKF
            big_R = kron(eye(length(yStoreNominal{ii})/3),Rtrue);
            K = phatm*H{ii}'*inv(H{ii}*phatm*H{ii}' + big_R);
            xP(ii,:) = xHatm(ii,:) + (K*errorY{ii})';
            P_predict = (eye(4) - K*H{ii})*phatm;
            err_x = (xTrue(ii,:) - xP(ii,:))';

            S_k = H{ii}(1:3,:)*phatm*H{ii}(1:3,:)' + Rtrue;                 % innovation covariance
            S_k = 0.5*(S_k + S_k');

            % NEES
            NEES(ii) = err_x'*(P_predict\err_x);

            % NIS
            NIS(ii) = errorY{ii}(1:3)'*(S_k\errorY{ii}(1:3));
        else
            % if we don't have measurements, consider xhatm as best estimate
            xP(ii,:) = xHatm(ii,:);
            P_predict = phatm;
        end

        stDev(ii,:) = 2*sqrt(diag(P_predict));
    end

    samplesNEES(n,:)    = NEES;
    samplesNIS(n,:)     = NIS;
    xTrueStore{n}       = xTrue;
end
%% Chi Square Testing


% NEES test
for zz = 2:length(tVec)
    nonzero = find(samplesNEES(:, zz) ~= 0);
    NEES_single = samplesNEES(:, zz);
    if length(nonzero) == simsN
        NEES_filtered(zz - 1) = mean(NEES_single(NEES_single ~= 0));
    end
end

NEES_filtered = NEES_filtered(find(NEES_filtered ~= 0));

alpha_NEES = 0.05;
Nn_x = simsN*length(F);

r1_x = chi2inv(alpha_NEES/2, Nn_x )./ simsN;
r2_x = chi2inv(1-alpha_NEES/2, Nn_x )./ simsN;

% NIS test
for zz = 2:length(tVec)
    nonzero2 = find(samplesNIS(:, zz) ~= 0);
    NIS_single = samplesNIS(:, zz);
    if length(nonzero2) == simsN
        NIS_filtered(zz - 1) = mean(NIS_single(NIS_single ~= 0));
    end
end

NIS_filtered = NIS_filtered(find(NIS_filtered ~= 0));

alpha_NIS = 0.05;
Nn_y = simsN*3;

r1_y = chi2inv(alpha_NIS/2, Nn_y )./ simsN;
r2_y = chi2inv(1-alpha_NIS/2, Nn_y )./ simsN;


tplot = linspace(0,steps*deltaT,steps);           % time plotting vector
NEES_time = linspace(1,length(tVec),length(tVec));      % time step plotting vector
NEES_time = NEES_time(find(NEES_filtered ~= 0));
NIS_time = linspace(1,length(tVec),length(tVec));
NIS_time = NIS_time(find(NIS_filtered ~= 0));


%% Plots
% NEES plot
figure
hold on
scatter(NEES_time,NEES_filtered)
plot(r1_x*ones(size(NEES_filtered)),'k--')
plot(r2_x*ones(size(NEES_filtered)),'k--')
title('NEES EKF')
xlabel('Time Step')
yl1 = ylabel('NEES statistic $\bar{\epsilon}_x$');
legend('NEES @ time k', 'r_1 bound', 'r_2 bound')
set(yl1,'interpreter','latex')

%
% % NIS plot
figure
hold on
scatter(NIS_time,NIS_filtered)
plot(r1_y*ones(size(NIS_filtered)),'k--')
plot(r2_y*ones(size(NIS_filtered)),'k--')
title('NIS EKF Results')
xlabel('Time Step')
yl2 = ylabel('NEES statistic $\bar{\epsilon}_y$');
legend('NIS @ time k', 'r_1 bound', 'r_2 bound')
set(yl2,'interpreter','latex')



%% Plots
% Noisy ground truth states
figure
hold on
grid on;
plot(tplot,xTrue(:,1))
plot(tplot,xTrue(:,3))
title('Noisy Ground Truth Position with Perturbation')
xlabel('Time (seconds)')
ylabel('State Positions (km)')
legend('X Position','Y Position')

figure
hold on
grid on;
plot(tplot,xTrue(:,2))
plot(tplot,xTrue(:,4))
title('Noisy Ground Truth Velocity with Perturbation')
xlabel('Time (seconds)')
ylabel('State Velocities (km/s)')
legend('X Velocity','Y Velocity')

figure
hold on
plot(tplot,xTrue(:,1)-xNom(:,1))
plot(tplot,xTrue(:,3)-xNom(:,3))
title('Difference of Noisy Ground Truth Position and Nominal State')
xlabel('Time (seconds)')
ylabel('State Positions (km)')
legend('X Position','Y Position')

figure
hold on
plot(tplot,xTrue(:,2)-xNom(:,2))
plot(tplot,xTrue(:,4)-xNom(:,4))
title('Difference of Noisy Ground Truth Velocity and Nominal State')
xlabel('Time (seconds)')
ylabel('State Velocities (km/s)')
legend('X Velocity','Y Velocity')

% Noisy measurements
[noisy_yplots] = plotYData(yStoreNoise, StationNumber, tVec);

% Linear KF state estimation errors
figure
hold on
title('Extended Kalman Filter Position Estimation Error')
plot(tplot,xP(:,1)'-xTrue(:,1)')
plot(tplot,xP(:,3)'-xTrue(:,3)')
xlabel('Time (seconds)')
ylabel('Position Error (km)')
legend('X Velocity','Y Velocity')

figure
hold on
title('Extended Kalman Filter Velocity Estimation Error')
plot(tplot,xP(:,2)'-xTrue(:,2)')
plot(tplot,xP(:,4)'-xTrue(:,4)')
xlabel('Time (seconds)')
ylabel('Velocity Error (km/s)')
legend('X Velocity','Y Velocity')



figure
subplot(4,1,1);
plot(tplot,xP(:,1));
hold on; grid on;
plot(tplot,xP(:,1) - stDev(1,:));
plot(tplot,xP(:,1) + stDev(1,:));
xlabel('Time (seconds)'); ylabel('X [km]');
legend('Trajectory','-2\sigma','2\sigma');
title('EKF Estimated States');

subplot(4,1,2);
plot(tplot,xP(:,2));
hold on; grid on;
plot(tplot,xP(:,2) - stDev(2,:));
plot(tplot,xP(:,2) + stDev(2,:));
xlabel('Time (seconds)'); ylabel('$\dot{X} [km/s]$','Interpreter','latex');

subplot(4,1,3);
plot(tplot,xP(:,3));
hold on; grid on;
plot(tplot,xP(:,3) - stDev(3,:));
plot(tplot,xP(:,3) + stDev(3,:));
xlabel('Time (seconds)'); ylabel('Y [km]');

subplot(4,1,4);
plot(tplot,xP(:,4));
hold on; grid on;
plot(tplot,xP(:,4) - stDev(4,:));
plot(tplot,xP(:,4) + stDev(4,:));
xlabel('Time (seconds)'); ylabel('$\dot{Y} [km/s]$','Interpreter','latex');





figure
subplot(4,1,1);
plot(tplot,xP(:,1)'-xTrue(:,1)')
hold on; grid on;
xlabel('Time (seconds)'); ylabel('X [km]');
title('EKF Estimated Error');

subplot(4,1,2);
plot(tplot,xP(:,2)'-xTrue(:,2)')
hold on; grid on;
xlabel('Time (seconds)'); ylabel('$\dot{X} [km/s]$','Interpreter','latex');

subplot(4,1,3);
plot(tplot,xP(:,3)'-xTrue(:,3)')
hold on; grid on;
xlabel('Time (seconds)'); ylabel('Y [km]');

subplot(4,1,4);
plot(tplot,xP(:,4)'-xTrue(:,4)')
hold on; grid on;
xlabel('Time (seconds)'); ylabel('$\dot{Y} [km/s]$','Interpreter','latex');


