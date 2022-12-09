%% full nonlinear sim
% constants
mu = 398600;
mu = 3.986004418e5;
Re = 6378;
we = 2*pi/86400;

% radius
r = @(x) sqrt(x(1,:).^2+x(3,:).^2);
% odefun
f = @(t,x) [x(2,:);...
           -mu*x(1,:)./r(x).^3;...
            x(4,:);...
           -mu*x(3,:)./r(x).^3];
% sim setup
dt = 10;
times = 0:dt:14000;
x_pert = [0;0.075;0;-0.021];
x_pert = [0;0;0;0];
x0 = [6678; 0; 0; sqrt(mu/6678)] + x_pert;

% solve ODE
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,x] = ode45(f,times,x0,opts);
t = t'; x = x'; i = 1:12; nt = length(t); ni = length(i);

% theta i
th = @(t,i) mod(we*t+(i-1)*pi/6+pi,2*pi) - pi;

% Xsi and d/dt Xsi
xi = @(t,i) Re*cos(we*t+(i-1)*pi/6); % t: row, i:col, xi: ni x nt
xi_dot = @(t,i) -we*Re*sin(we*t+(i-1)*pi/6);
% Ysi and d/dt Ysi
yi = @(t,i) Re*sin(we*t+(i-1)*pi/6); % t: row, i:col, yi: ni x nt
yi_dot = @(t,i) we*Re*cos(we*t+(i-1)*pi/6);

% range measurement
p = @(t,x,i) sqrt((x(1,:)-xi(t,i)).^2 + (x(3,:)-yi(t,i)).^2); % t: row, x: 4 x nt, i: col, p: ni x nt
% range rate measurement
p_dot = @(t,x,i) ( (x(1,:)-xi(t,i)).*(x(2,:)-xi_dot(t,i)) + (x(3,:)-yi(t,i)).*(x(4,:)-yi_dot(t,i)) )./p(t,x,i);
% angle measurement
phi = @(t,x,i) atan2(x(3,:)-yi(t,i), x(1,:)-xi(t,i));
% total measurement
h = @(t,x,i) [p(t,x,i);...
              p_dot(t,x,i);...
              phi(t,x,i)];

% calculate measurements for all 12 ground stations
nonlin_y = zeros(3,nt,ni);
figure(33); clf(); colors = 'rgbcmkrgbcmk';
for i = 1:ni
    % angle measurements over time
    phis = phi(t,x,i);
    % ground station position over time
    ths = th(t,i);

    % vector measurement over time for ground station i
    nonlin_y(:,:,i) = [p(t,x,i); p_dot(t,x,i); phis];

    % not in view when phi is more than 90 deg from theta
    not_in_view = min(mod(ths-phis,2*pi),mod(phis-ths,2*pi)) > pi/2;

    % set measurements to NaN when not in view
    nonlin_y(:,not_in_view,i) = NaN;

    % plot measurements
    subplot(3,1,1); plot(t,nonlin_y(1,:,i),[colors(i) 'o']); hold on
    subplot(3,1,2); plot(t,nonlin_y(2,:,i),[colors(i) 'o']); hold on
    subplot(3,1,3); plot(t,nonlin_y(3,:,i),[colors(i) 'o']); hold on
end

%% linearization

% CT jacobians
A = @(x) [0,                             1,                              0, 0; ...
         -mu/r(x)^3*(1-3*x(1)^2/r(x)^2), 0,          3*mu*x(1)*x(3)/r(x)^5, 0;...
          0,                             0,                              0, 1;...
          3*mu*x(1)*x(3)/r(x)^5,         0, -mu/r(x)^3*(1-3*x(3)^2/r(x)^2), 0];
B = [0 0;1 0;0 0;0 1];

H = @(t,x,i) [(x(1)-xi(t,i))/p(t,x,i), 0, (x(3)-yi(t,i))/p(t,x,i), 0;...
              (x(2)-xi_dot(t,i))/p(t,x,i) - p_dot(t,x,i)*(x(1)-xi(t,i))/p(t,x,i)^2,...
                (x(1)-xi(t,i))/p(t,x,i),...
                (x(4)-yi_dot(t,i))/p(t,x,i) - p_dot(t,x,i)*(x(3)-yi(t,i))/p(t,x,i)^2,...
                (x(3)-yi(t,i))/p(t,x,i);...
              cos(phi(t,x,i))^2*[-(x(3)-yi(t,i))/(x(1)-xi(t,i))^2, 0, 1/(x(1)-xi(t,i)), 0]];

% set up nominal trajectory
% nominal radius/mean motion/period
r0 = 6678; n0 = sqrt(mu/r0^3); T0 = 2*pi/n0;
x_nom = @(t) [r0*cos(n0*t);...
            -r0*n0*sin(n0*t);...
             r0*sin(n0*t);...
             r0*n0*cos(n0*t)];

% simulate using DT jacobian
Alin_CT = zeros(4,4,nt);
Alin_DT = Alin_CT;
x_lin = zeros(4,nt);
dx = zeros(4,nt+1);
dx(:,1) = x_pert;
y_lin = zeros(3,nt,ni);
dy = zeros(3,nt,ni);
for k = 1:nt
    % current state of nominal trajectory
    x_nom_k = x_nom(t(k));

    % reconstruct full state as nominal + delta
    x_lin(:,k) = x_nom_k + dx(:,k);
    
    % calculate linearized measurement
    for i = 1:ni
        % evaluate nominal measurement
        y_lin(:,k,i) = [p(t(k),x_nom_k,i); p_dot(t(k),x_nom_k,i); phi(t(k),x_nom_k,i)];

        % add delta measurement
        y_lin(:,k,i) = y_lin(:,k,i) + H(t(k),x_nom_k,i)*dx(:,k);

        % remove if out of sight
        if min(mod(th(t(k),i)-phi(t(k),x_lin(:,k),i),2*pi),...
               mod(phi(t(k),x_lin(:,k),i)-th(t(k),i),2*pi)) > pi/2
            y_lin(:,k,i) = NaN;
        end
    end

    % evaluate CT jacobian at the nominal trajectory
    Alin_CT(:,:,k) = A(x_nom_k);
    
    % euler estimate of DT jacobian
    Alin_DT(:,:,k) = eye(4) + dt*Alin_CT(:,:,k);

    % propagate the delta x state
    dx(:,k+1) = Alin_DT(:,:,k)*dx(:,k);
end
G = dt*B;
W = G;
I = eye(4);


% plot linearized sim state on top of the full nonlinear state
figure(44); clf()
subplot(4,1,1); plot(t,x(1,:));
hold on; plot(t,x_lin(1,:));
subplot(4,1,2); plot(t,x(2,:));
hold on; plot(t,x_lin(2,:));
subplot(4,1,3); plot(t,x(3,:));
hold on; plot(t,x_lin(3,:));
subplot(4,1,4); plot(t,x(4,:));
hold on; plot(t,x_lin(4,:));

% plot the difference of the two sims
figure(55); clf()
subplot(4,1,1); plot(t,x(1,:)-x_lin(1,:));
subplot(4,1,2); plot(t,x(2,:)-x_lin(2,:));
subplot(4,1,3); plot(t,x(3,:)-x_lin(3,:));
subplot(4,1,4); plot(t,x(4,:)-x_lin(4,:));

% plot linearized sim state on top of the full nonlinear state
figure(66); clf(); dx = dx(:,1:end-1);
subplot(4,1,1); plot(t,dx(1,:));
subplot(4,1,2); plot(t,dx(2,:));
subplot(4,1,3); plot(t,dx(3,:));
subplot(4,1,4); plot(t,dx(4,:));

% pot linearized measurements
figure(77); clf()
colors = 'rgbcmkrgbcmk';
for i = 1:ni
    subplot(3,1,1); plot(t,y_lin(1,:,i),[colors(i) 'o']); hold on
    subplot(3,1,2); plot(t,y_lin(2,:,i),[colors(i) 'o']); hold on
    subplot(3,1,3); plot(t,y_lin(3,:,i),[colors(i) 'o']); hold on
end

%% generate noisy measurements


%% Linearized KF
load('orbitdeterm_finalproj_KFdata.mat');

% noise parameters
Q = Qtrue;
R = 10*Rtrue;

x0 = [0;0;0;0];
P0 = 300 * diag([1,1e-3,1,1e-3]);
[xs_lkf,Ps_lkf,s_lkf,inns_lkf] = LKF(x0,P0,x_nom,h,A,H,W,Q,R,t,ydata);

%% EKF

x0 = x_nom(0);
P0 = 300 * diag([1,1e-3,1,1e-3]);
h = @(t,x,i) [p(t,x,i);...
              p_dot(t,x,i);...
              phi(t,x,i)];
[xs_ekf,Ps_ekf,s_ekf,inns_ekf] = EKF(x0,P0,f,h,A,H,W,Qtrue,Rtrue,t,ydata);

%% filter functions

function [xs,Ps,s,inns] = LKF(x0,P0,x_nom,h,A,H,W,Q,R,t,ydata)
% initial conditions
dt = t(2)-t(1);
ny = length(ydata);
x = x0;
P = P0;
xs = zeros(4,ny); xs(:,1) = x;
Ps = zeros(4,4,ny); Ps(:,:,1) = P;
s = zeros(4,ny); s(:,1) = 2*sqrt(diag(P));
inns = zeros(3,ny);
I = eye(length(x0));
for k = 2:ny
    % get model parameters
    x_nom_k = x_nom(t(k));
    Fk = I + dt*A(x_nom(t(k-1)));
    if ~isempty(ydata{k})
        i_in_view = ydata{k}(4,:);
        yk = reshape(ydata{k}(1:3,:),[],1);
    else
        i_in_view = [];
        yk = [];
    end
    Hs = cell(1,length(i_in_view)); Rs = Hs;
    y_nom_k = zeros(3*length(i_in_view),1);
    for j = 1:length(i_in_view)
        Hs{j} = H(t(k),x_nom_k,i_in_view(j));
        Rs{j} = R;
        y_nom_k((j-1)*3+(1:3)) = h(t(k),x_nom_k,i_in_view(j));
    end
    Hk = vertcat(Hs{:});
    Rk = blkdiag(Rs{:});
    
    % predict
    x_ = Fk*x;
    P_ = Fk*P*Fk' + W*Q*W';
    % correct
    if ~isempty(i_in_view)
        K = P_*Hk'/(Hk*P_*Hk' + Rk);
        innovation = (yk-y_nom_k) - Hk*x_;
        innovation(3:3:end) = mod(pi + innovation(3:3:end), 2*pi) - pi;
        x = x_ + K*innovation;
        P = (I-K*Hk)*P_;
        inns(:,k) = innovation(1:min(end,3));
    else
        x = x_;
        P = P_;
    end
    
    % store current values
    xs(:,k) = x;
    Ps(:,:,k) = P;
    s(:,k) = 2*sqrt(diag(P));
end
% reconstruct full state estimate 
xs = x_nom(t) + xs;
end



function [xs,Ps,s,inns] = EKF(x0,P0,f,h,A,H,W,Q,R,t,ydata)
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
ny = length(ydata);
dt = t(2)-t(1);
x = x0;
P = P0;
I = eye(length(x0));
xs = zeros(4,ny); xs(:,1) = x;
Ps = zeros(4,4,ny); Ps(:,:,1) = P;
s = zeros(4,ny); s(:,1) = 2*sqrt(diag(P));
inns = zeros(3,ny);
for k = 2:ny
    % get F
    Fk = I + dt*A(x);
    
    % prediction step
    [~,x_] = ode45(f,[0 dt],x,opts);
    x_ = x_(end,:)';
    P_ = Fk*P*Fk' + W*Q*W';

    % measurement processing
    if ~isempty(ydata{k})
        i_in_view = ydata{k}(4,:);
        yk = reshape(ydata{k}(1:3,:),[],1);
    else
        i_in_view = [];
        yk = [];
    end
    Hs = cell(1,length(i_in_view)); Rs = Hs;
    yk_est = zeros(3*length(i_in_view),1);
    for j = 1:length(i_in_view)
        Hs{j} = H(t(k),x,i_in_view(j));
        Rs{j} = R;
        yk_est((j-1)*3+(1:3)) = h(t(k),x_,i_in_view(j));
    end

    % create augmented H and R matrices
    Hk = vertcat(Hs{:});
    Rk = blkdiag(Rs{:});

    % correction step
    if ~isempty(i_in_view)
        K = P_*Hk'/(Hk*P_*Hk' + Rk);
        innovation = yk - yk_est;
        innovation(3:3:end) = mod(pi + innovation(3:3:end), 2*pi) - pi;
        x = x_ + K*innovation;
        P = (I-K*Hk)*P_;
        inns(:,k) = innovation(1:min(end,3));
    else
        x = x_;
        P = P_;
    end

    % store current values
    xs(:,k) = x;
    Ps(:,:,k) = P;
    s(:,k) = 2*sqrt(diag(P));
end
end

