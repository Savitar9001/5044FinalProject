%% full nonlinear sim
% constants
mu = 398600;
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

% calculate measurements for all 12 ground stations
h = zeros(3,nt,ni);
figure(33); clf()
for i = 1:ni
    % angle measurements over time
    phis = phi(t,x,i);
    % ground station position over time
    ths = th(t,i);

    % vector measurement over time for ground station i
    h(:,:,i) = [p(t,x,i); p_dot(t,x,i); phis];

    % not in view when phi is more than 90 deg from theta
    not_in_view = min(mod(ths-phis,2*pi),mod(phis-ths,2*pi)) > pi/2;

    % set measurements to NaN when not in view
    h(:,not_in_view,i) = NaN;

    % plot measurements
    subplot(3,1,1); plot(t,h(1,:,i),'o'); hold on
    subplot(3,1,2); plot(t,h(2,:,i),'o'); hold on
    subplot(3,1,3); plot(t,h(3,:,i),'o'); hold on
end

%% linearization

% CT jacobians
A = @(x) [0, 1, 0, 0; ...
         -mu/r(x)^3*(1+3*x(1)^2/r(x)^2), 0, 3*mu*x(1)*x(3)/r(x)^5, 0;...
          0, 0, 1, 0;...
          3*mu*x(1)*x(3)/r(x)^5, 0, -mu/r(x)^3*(1+3*x(3)^2/r(x)^2), 0];
B = [0 0;1 0;0 0;0 1];

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
for k = 1:nt
    % reconstruct full state as nominal + delta
    x_lin(:,k) = x_nom(t(k)) + dx(:,k);

    % evaluate CT jacobian at the nominal trajectory
    Alin_CT(:,:,k) = A(x_nom(t(k)));
    
    % euler estimate of DT jacobian
    Alin_DT(:,:,k) = eye(4) + dt*Alin_CT(:,:,k);

    % propagate the delta x state
    dx(:,k+1) = Alin_DT(:,:,k)*dx(:,k);
end

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

% i = 2;
% % h: ny x nt x ni
% h = @(t,x,i) [reshape(p(t,x,i)',1,nt,ni);...
%               reshape(p_dot(t,x,i)',1,nt,ni);...
%               reshape(phi(t,x,i)',1,nt,ni)];