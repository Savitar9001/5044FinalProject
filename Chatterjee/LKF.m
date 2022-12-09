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
state0              = 6678;                                 % [km] Nominal Orbit of SC
y0              = 0;
mu              = 398600;                               % [km^3/s^2] Earth Gravitational Parameters
r0              = sqrt(state0^2 + y0^2);
yDot0           = r0*sqrt(mu/(r0^3));
rE              = 6378;                                 % [km] Earth Radius Equatorial
omega           = sqrt(mu/r0^3);
spinRateEarth   = (2*pi)/86400;                         % [rad/s] Spinrate of the Earth
deltaT          = 10;
dtVec           = [0,10];
Per             = ceil(2*pi*sqrt(r0^3/mu));
steps           = 1401;
tVec            = linspace(0,14000,steps);              % Time Vector
TA              = linspace(0,2*pi,steps);               % [rad] True Anamoly
stationNum      = 12;
trans           = [0 1 0;-1 0 0; 0 0 1];
simsN           = 50;

% Atest           = A_nom([6678, 0, 0, yDot0]); % Use for Debugging


%% Start Prelminary Calculations
theta0 = NaN(stationNum,1);     % Preallocate to Speed Up
for ii = 1:stationNum
    theta0(ii) = (ii-1)*(pi/6);                      % Calculate Distribution of the Stations
end

gamma = [0 0; 1 0; 0 0; 0 1];                       % Disturbance model
tuneQ = diag([1e-8,1e-3,1e-8,1e-3]);               % just chose small #s for now
%Q_tune = diag([10,0.1,10,0.1]);
%Omega_ext = deltaT*Gamma;

perturb = [0.01, 0.075, 0.01, -0.021];

% Nonlinear Simulation
state0  = [r0*cos(omega*t0), -omega*r0*sin(omega*t0), r0*sin(omega*t0), omega*r0*cos(omega*t0)] + perturb;
P0      = diag((perturb.^2));                 % Might Need to Double Check This
P       = 2*pi*sqrt(r0^3/mu);                 % [sec] Period

tol = 1e-12;
options = odeset('RelTol',tol,'AbsTol',[tol tol tol tol]);

[nonLinearT,nonLinearState] = ode45(@propDyDt,tVec,state0,options);


%Rtrue = 0.68*Rtrue;

%% The Following Are Trials of the Simulation
samplesNEES = zeros(simsN,length(tVec));
samplesNIS = zeros(simsN,length(tVec));
for n = 1:simsN

    trueX = nonLinearState(1,:);                          % Initial estimate

    % simulation with process noise
    for ii = 2:length(tVec)
        wk = mvnrnd(zeros(1,2),Qtrue); % Qtrue is loaded from the mat file
        [nonLinearT,noisyX] = ode45(@propDyDTNoise,dtVec,trueX(ii-1,:),options);
        trueX(ii,:) = noisyX(end,:);
    end

    % Let's Iterate Through Time
    H               = cell(length(tVec)-1,1);                   % Let's Try to Preallocate to be a Good Programmer
    HExt            = cell(length(tVec)-1,1);
    storeF          = cell(length(tVec)-1,1);
    yStackNom       = NaN(length(tVec), 3*stationNum);
    yStackNoise     = NaN(length(tVec), 3*stationNum);
    StationNumber   = NaN(stationNum,length(tVec));
    yStoreNominal   = cell(length(tVec),1);
    yStoreNoise     = cell(length(tVec),1);
    xStation    = NaN(length(tVec), 2*stationNum);
    yStation    = NaN(length(tVec), 2*stationNum);

    for ii = 2:length(tVec)
        A = findANominal(nonLinearState(ii-1,:));                   % CT system matrix at each time step
        F = eye(4) + deltaT*A;                                      % DT system matrix at each time step
        storeF{ii-1} = F;

        %         % EKF Matrices
        %         A_ext = A_nom(trueX(ii-1,:));
        %         F_ext = eye(4) + deltaT*A_ext;
        %         F_ext_store{ii-1} = F;

        tempX       = NaN(1,2*stationNum);
        tempY       = NaN(1,2*stationNum);
        temp_yOut   = cell(stationNum,1);
        temp_yNoise = cell(stationNum,1);
        % yNominal = NaN(length(tVec), 3*stationNum);
        H{ii}       = [];
        HExt{ii}    = [];

        for kk = 1:stationNum
            tempX(2*kk-1:2*kk)  = [rE*cos(spinRateEarth*ii*deltaT + theta0(kk));-spinRateEarth*rE*sin(spinRateEarth*ii*deltaT+theta0(kk))];
            tempY(2*kk-1:2*kk)  = [rE*sin(spinRateEarth*ii*deltaT + theta0(kk));spinRateEarth*rE*cos(spinRateEarth*ii*deltaT+theta0(kk))];

            temp_yOut{kk,1}     = [measureY(nonLinearState(ii,:),tempX(2*kk-1:2*kk),tempY(2*kk-1:2*kk))];
            temp_yNoise{kk,1}   = measureY(trueX(ii,:),tempX(2*kk-1:2*kk),tempY(2*kk-1:2*kk))+mvnrnd(zeros(1,3),Rtrue);
        end
        temp_yOut   = [temp_yOut{:}]; % Converting Cell to Matrix
        temp_yNoise = [temp_yNoise{:}];

        xStation(ii,:)      = tempX;
        yStation(ii,:)      = tempY;
        yNominal            = temp_yOut(:);
        y_noise             = temp_yNoise(:);
        yStackNom(ii,:)     = yNominal;
        yStackNoise(ii,:)   = y_noise;

        tempStation = zeros(12,1);

        for jj = 1:stationNum
            pos = [xStation(ii,2*jj-1),yStation(ii,2*jj-1),0];
            pos2 = [trueX(ii,1) - xStation(ii,2*jj-1),trueX(ii,3)-yStation(ii,2*jj-1),0];
            posh = trans*pos';
            c = cross(posh,pos2);

            if c(3) >= 0
                tempStation(jj) = jj;
                H{ii} = [H{ii};findHNominal(nonLinearState(ii,:),tempX(2*jj-1:2*jj),tempY(2*jj-1:2*jj))];

                % EKF matrix
                HExt{ii} = [HExt{ii};findHNominal(trueX(ii,:),tempX(2*jj-1:2*jj),tempY(2*jj-1:2*jj))];
            end

            if c(3) < 0
                yNominal(3*jj-2:3*jj) = 0;
                y_noise(3*jj-2:3*jj) = 0;
            end

            yNominal(3*jj) = wrapToPi(yNominal(3*jj));
            y_noise(3*jj) = wrapToPi(y_noise(3*jj));

        end

        StationNumber(:,ii) = tempStation;
        yStoreNominal{ii} = yNominal(find(yNominal~=0));
        yStoreNoise{ii} = y_noise(find(y_noise~=0));
    end

    %    kfdx = []; P = []; phatm = []; dy = []; dy_pert = []; stdev = []; innov = [];

    [kfdx, P, phatm, dy, dy_pert, stdev, innov] = linearKalmanFilter(yStoreNoise, yStoreNominal, perturb, P0, storeF, gamma, tuneQ, Rtrue, H, deltaT);

    %% NEES and NIS Testing

    % Preallocate Storage
    errorY      = cell(length(tVec)-1,1);
    NEES        = zeros(length(tVec),1);
    NIS         = zeros(length(tVec),1);

    for kk = 2:length(tVec)
        if length(yStoreNoise{kk}) > 1 % make sure we have data

            errorY{kk} = (yStoreNoise{kk} - yStoreNominal{kk}) - dy_pert{kk};

            for z = 1:length(errorY{kk})/3 % wrap per measurement
                tempErrorY = errorY{kk};
                tempErrorY(3*z) = wrapToPi(tempErrorY(3*z));
                errorY{kk} = tempErrorY;
            end

            errorX = (trueX(kk,:) - (nonLinearState(kk,:) - kfdx(kk,:)))';

            P2 = 0.5*(P(:,:,kk) + P(:,:,kk)');
            SK = H{kk}(1:3,:)*phatm(:,:,kk)*H{kk}(1:3,:)' + Rtrue;                 % innovation covariance
            SK = 0.5*(SK + SK');

            % NEES
            NEES(kk) = errorX'*(P2\errorX);

            % NIS
            NIS(kk) = errorY{kk}(1:3)'*(SK\errorY{kk}(1:3));
        end
    end

    samplesNEES(n,:)    = NEES;
    samplesNIS(n,:)     = NIS;
    NEES                = zeros(length(tVec),1);
    NIS                 = zeros(length(tVec),1);
end


H3      = cell(1,length(tVec));
yNom3   = cell(1,length(tVec));
y3      = cell(1,length(tVec));
for ii = 2:length(tvec)
    if size(ydata{ii},2) == 2
        y3{ii} = [ydata{ii}(1:3,1); ydata{ii}(1:3,2)];
        station1 = ydata{ii}(4,1);
        station2 = ydata{ii}(4,2);
        yNom3{ii} = [yStackNom(ii,station1*3-2:station1*3)'; yStackNom(ii,station2*3-2:station2*3)'];
        H3{ii} = [findHNominal(nonLinearState(ii,:),xStation(ii,2*station1-1:2*station1),yStation(ii,2*station1-1:2*station1));
            findHNominal(nonLinearState(ii,:),xStation(ii,2*station2-1:2*station2),yStation(ii,2*station2-1:2*station2))];
    elseif size(ydata{ii},2) == 1
        y3{ii} = ydata{ii}(1:3);
        station1 = ydata{ii}(4);
        yNom3{ii} = yStackNom(ii,station1*3-2:station1*3)';
        H3{ii} = [findHNominal(nonLinearState(ii,:),xStation(ii,2*station1-1:2*station1),yStation(ii,2*station1-1:2*station1))];
    else
        y3{ii} = ydata{ii};
        yNom3{ii} = [];
        H3{ii} = [];
    end
end

% kfdx_p3 = []; P_p3 = []; phatm_p3 = []; dy_p3 = []; dy_pert_p3 = []; stdev_p3 = []; innov_p3 = [];
[kfdx_p3, P_p3, phatm_p3, dy_p3, dy_pert_p3, stdev_p3] = linearKalmanFilter(y3, yNom3, zeros(1,4), P0, storeF, gamma, tuneQ, Rtrue, H3, deltaT);



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
title('NEES LKF')
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
title('NIS LKF Results')
xlabel('Time Step')
yl2 = ylabel('NEES statistic $\bar{\epsilon}_y$');
legend('NIS @ time k', 'r_1 bound', 'r_2 bound')
set(yl2,'interpreter','latex')


%% Plots
% Noisy ground truth states
figure
hold on
grid on;
plot(tplot,trueX(:,1))
plot(tplot,trueX(:,3))
title('Noisy Ground Truth Position with Perturbation')
xlabel('Time [sec]')
ylabel('State Positions [km]')
legend('X','Y')

figure
hold on
grid on
plot(tplot,trueX(:,2))
plot(tplot,trueX(:,4))
title('Noisy Ground Truth Velocity with Perturbation')
xlabel('Time [sec]')
ylabel('State Velocities [km/s]')
legend('X','Y')

figure
hold on
grid on
plot(tplot,trueX(:,1)-nonLinearState(:,1))
plot(tplot,trueX(:,3)-nonLinearState(:,3))
title('Difference of Noisy Ground Truth Position and Nominal State')
xlabel('Time [sec]')
ylabel('State Positions [km]')
legend('X','Y')

figure
hold on
grid on
plot(tplot,trueX(:,2)-nonLinearState(:,2))
plot(tplot,trueX(:,4)-nonLinearState(:,4))
title('Difference of Noisy Ground Truth Velocity and Nominal State')
xlabel('Time [sec]')
ylabel('State Velocities [km/s]')
legend('X','Y')

% Noisy measurements
[noisy_yplots] = plotYData(yStoreNoise, StationNumber, tVec);

% Linear KF state estimation errors
figure
hold on
grid on
title('LKF Position Est Error')
plot(tplot,kfdx(:,1)')
plot(tplot,kfdx(:,3)')
xlabel('Time (seconds)')
ylabel('Position Error [km]')
legend('X','Y')

figure
hold on
grid on
title('LFK Velocity Est Error')
plot(tplot,kfdx(:,2)')
plot(tplot,kfdx(:,4)')
xlabel('Time [sec]')
ylabel('Velocity Error [km/s]')
legend('X','Y')

