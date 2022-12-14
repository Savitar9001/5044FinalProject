clearvars -except xs y
% [xs_ekf6,Ps_ekf6,s_ekf6,invS6,dy6,inns_ekf6] = EKF(x0,P0,f,h,A,H,W,Q_ekf,Rtrue,t,ydata);
%[xs_sim,y_sim] = sim([6704.9;0.1;-112.9;7.8],t,Gamma,Q,R,f,p,p_dot,phi,th);

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
%x_pert = [0;0;0;0];
x0 = [6678; 0; 0; sqrt(mu/6678)] + x_pert;

% solve ODE
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,x] = ode45(f,times,x0,opts);
t = t'; x = x'; is = 1:12; nt = length(t); ni = length(is);

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
Gamma = B;
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

% % Floquent analysis of A stability
% dPhi_dt = @(t,x) reshape( A(x_nom(t))*reshape(x,4,4), 16,1);
% [~,Phis] = ode45(dPhi_dt,[0,T0],reshape(eye(4),16,1),opts);
% Phi = reshape(Phis(end,:)',4,4);
% dPhi_dt = @(x) reshape( A(x(1:4))*reshape(x(5:20),4,4), 16,1);
% dXdt = @(t,x)
% [~,Phis] = ode45(dPhi_dt,[0,T0],reshape(eye(4),16,1),opts);
% Phi = reshape(Phis(end,:)',4,4);

%dXdt = @(t,x) [f(t,x(1:4)); dPhidt(t,x)];

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
P0 = 50^2 * diag([1,1e-6,1,1e-6]);
P0 = diag([0.0001 0.0056 0.0001 0.0004]);
%% generate noisy measurements
load('orbitdeterm_finalproj_KFdata.mat');
Q = Qtrue;
R = Rtrue;
% x_pert = chol(P0,'lower')*randn(4,1) 
% [xs_sim,y_sim] = sim(x_nom(0)+x_pert,t,Gamma,Q,R,f,p,p_dot,phi,th);
% 
% %% Linearized KF
% 
% x0 = [0;0;0;0];
% [xs_lkf,Ps_lkf,s_lkf,inns_lkf] = LKF(x0,P0,x_nom,h,A,H,W,Q,R,t,y_sim);
% 
% %% EKF
% 
% x0 = x_nom(0);
% %P0 = 3000 * diag([1,1e-3,1,1e-3]);
% [xs_ekf,Ps_ekf,s_ekf,inns_ekf] = EKF(x0,P0,f,h,A,H,W,Qtrue,Rtrue,t,y_sim);


%% NEES/NIS tests
num_runs = 100;
regen_TMT_data = 0;
alpha = 0.05;
if regen_TMT_data
    clear xs y
    xs = cell(1,num_runs);
    y = cell(1,num_runs);
    for j = 1:num_runs
        j
        x0 = x_nom(0) + chol(P0,'lower')*randn(4,1)
        [xs{j},y{j}] = sim(x0,t,Gamma,Q,R,f,p,p_dot,phi,th);
    end
end

% LKF
Q_lkf = 1*Q;
Q_lkf = 0.01*[0.0071246231041958     0.001866095612661221e-2;
         0.001866095612661221e-2    0.00779388288435167];
W(1,1) = 0.5*dt^2; W(3,2) = 05*dt^2; 
Q_lkf= W*Qtrue*W';
x0 = [0;0;0;0];
%P0 = 50^2 * diag([1,1e-6,1,1e-6]);
for j = 1:num_runs
    [xs_lkf(:,:,j),Ps_lkf(:,:,:,j),s_lkf(:,:,j),invS,dy,~] = LKF(x0,P0,x_nom,h,A,H,Q_lkf,Rtrue,t,y{j});
    ex(j,:) = pagemtimes(pagemtimes(permute(xs{j}-xs_lkf(:,:,j),[3 1 2]),pageinv(Ps_lkf(:,:,:,j))),...
                         permute(xs{j}-xs_lkf(:,:,j),[1 3 2]));
    for k = 1:nt
        numel_y(j,k) = sum(~isnan(dy{k}));
        ey(j,k) = dy{k}'*invS{k}*dy{k};
    end
    
end
exn = sum(ex)./num_runs;
r1 = chi2inv(alpha/2,num_runs*4)/num_runs;
r2 = chi2inv(1-alpha/2,num_runs*4)/num_runs;
figure(101);clf(); plot(t,exn,'.'); yline(r1,'--r'); yline(r2,'--r') 
title('LKF NEES'); axis([0 14000 0 2*r2(end)])
xlabel('time (seconds)'); ylabel('NEES Statistic')

py = sum(numel_y)/num_runs;
py = numel_y(end,:);
eyn = sum(ey.*~isnan(ey))./sum(~isnan(ey));
r1 = chi2inv(alpha/2,num_runs*py)./sum(~isnan(ey));%num_runs;
r2 = chi2inv(1-alpha/2,num_runs*py)./sum(~isnan(ey));%num_runs;
figure(102);clf(); plot(t,eyn,'.');hold on; plot(t,r1,'--r'); plot(t,r2,'--r');
title('LKF NIS'); axis([0 14000 0 2*r2(end)])
xlabel('time (seconds)'); ylabel('NIS Statistic')

%% EKF
Q_ekf = 1*Q;
x0 = x_nom(0);
W(1,1) = 0.5*dt^2; W(3,2) = 05*dt^2; 
Q_ekf = 1.4*W*Qtrue*W';
%P0 = 50^2 * diag([1,1e-3,1,1e-3]);
for j = 1:num_runs
    [xs_ekf(:,:,j),Ps_ekf(:,:,:,j),s_ekf(:,:,j),invS,dy,inns_ekf] = EKF(x0,P0,f,h,A,H,Q_ekf,Rtrue,t,y{j});
    ex(j,:) = pagemtimes(pagemtimes(permute(xs{j}-xs_ekf(:,:,j),[3 1 2]),pageinv(Ps_ekf(:,:,:,j))),...
                         permute(xs{j}-xs_ekf(:,:,j),[1 3 2]));
    for k = 1:nt
        numel_y(j,k) = length(dy{k});
        ey(j,k) = dy{k}'*invS{k}*dy{k};
    end
    
end
exn = sum(ex)./num_runs;
r1 = chi2inv(alpha/2,num_runs*4)/num_runs;
r2 = chi2inv(1-alpha/2,num_runs*4)/num_runs;
figure(103);clf(); plot(t,exn,'.'); yline(r1,'--r'); yline(r2,'--r');
title('EKF NEES');axis([0 14000 0 2*r2(end)])
xlabel('time (seconds)'); ylabel('NEES Statistic')

% eyn = sum(ey)./num_runs;
% r1 = chi2inv(alpha/2,num_runs*numel_y(1,:))/num_runs;
% r2 = chi2inv(1-alpha/2,num_runs*numel_y(1,:))/num_runs;
py = sum(numel_y)./sum(~isnan(ey));
%py = numel_y(end,:);
eyn = sum(ey.*~isnan(ey))./sum(~isnan(ey));
r1 = chi2inv(alpha/2,num_runs*py)./sum(~isnan(ey));%num_runs;
r2 = chi2inv(1-alpha/2,num_runs*py)./sum(~isnan(ey));%num_runs;
figure(104);clf(); plot(t,eyn,'.');hold on; plot(t,r1,'--r'); plot(t,r2,'--r');
title('EKF NIS'); axis([0 14000 0 2*r2(end)])
xlabel('time (seconds)'); ylabel('NIS Statistic')

%% part 6
[xs_ekf6,Ps_ekf6,s_ekf6,invS6,dy6,inns_ekf6] = EKF(x0,P0,f,h,A,H,Q_ekf,Rtrue,t,ydata);

%% functions

function [xs,y] = sim(x0,t,Gamma,Q,R,f,p,p_dot,phi,th)
opts = odeset('RelTol',1e-12,'AbsTol',1e-12); colors = 'rgbcmkrgbcmk';
dt = t(2)-t(1);
nt = length(t); is = 1:12; ni = length(is);
sqrtQ = chol(Q,'lower');
% simulate noisy dynamics
x = x0;
for j = 1:length(t)
    xs(:,j) = x;
    w = sqrtQ*randn(2,1);
    f_stoch = @(t,x) f(t,x) + Gamma*w;
    for k = 1:10
        [~,y] = ode45(f_stoch,[0,dt/10],x,opts);
        x = y(end,:)';
    end
end

% noisy measurements
ys = zeros(3,nt,ni);
i_in_view = false(ni,nt);
sqrtR = chol(R,'lower');
for i = 1:ni
    % angle measurements over time
    phis = phi(t,xs,i);
    % ground station position over time
    ths = th(t,i);
    
    % measurement noise v
    v = sqrtR*randn(3,nt);

    % vector measurement over time for ground station i
    ys(:,:,i) = [p(t,xs,i); p_dot(t,xs,i); phis] + v;

    % not in view when phi is more than 90 deg from theta
    not_in_view = min(mod(ths-phis,2*pi),mod(phis-ths,2*pi)) > pi/2;
    
    i_in_view(i,:) = min(mod(ths-phis,2*pi),mod(phis-ths,2*pi)) <= pi/2;
    
    % set measurements to NaN when not in view
    ys(:,not_in_view,i) = NaN;

    % plot measurements
    subplot(3,1,1); plot(t,ys(1,:,i),[colors(i) 'o']); hold on
    subplot(3,1,2); plot(t,ys(2,:,i),[colors(i) 'o']); hold on
    subplot(3,1,3); plot(t,ys(3,:,i),[colors(i) 'o']); hold on
end

y = cell(1,nt);
for j = 1:nt
    yj = reshape(squeeze(ys(:,j,:)),[],1);
    yj(isnan(yj)) = [];
    yj = reshape(yj,3,[]);
    ij = is(i_in_view(:,j));
    y{j} = [yj;ij];
end
end
 



% function [xs,Ps,s,invSs,dy,inns] = LKF(x0,P0,x_nom,h,A,H,W,Q,R,t,ydata)
% % initial conditions
% dt = t(2)-t(1);
% ny = length(ydata);
% x = x0;
% P = P0;
% xs = zeros(4,ny); xs(:,1) = x;
% Ps = zeros(4,4,ny); Ps(:,:,1) = P;
% invSs = cell(1,ny);
% dy = cell(1,ny);
% s = zeros(4,ny); s(:,1) = 2*sqrt(diag(P));
% inns = zeros(3,ny);
% I = eye(length(x0));
% for k = 2:ny
%     % get model parameters
%     x_nom_k = x_nom(t(k));
%     Fk = I + dt*A(x_nom(t(k-1)));
%     if ~isempty(ydata{k})
%         i_in_view = ydata{k}(4,:);
%         yk = reshape(ydata{k}(1:3,:),[],1);
%     else
%         i_in_view = [];
%         yk = [];
%     end
%     Hs = cell(1,length(i_in_view)); Rs = Hs;
%     y_nom_k = zeros(3*length(i_in_view),1);
%     for j = 1:length(i_in_view)
%         Hs{j} = H(t(k),x_nom_k,i_in_view(j));
%         Rs{j} = R;
%         y_nom_k((j-1)*3+(1:3)) = h(t(k),x_nom_k,i_in_view(j));
%     end
%     Hk = vertcat(Hs{:});
%     Rk = blkdiag(Rs{:});
%     
%     % predict
%     x_ = Fk*x;
%     P_ = Fk*P*Fk' + W*Q*W';
%     % correct
%     if ~isempty(i_in_view)
%         invS = inv(Hk*P_*Hk' + Rk);
%         K = P_*Hk'/(Hk*P_*Hk' + Rk);
%         innovation = (yk-y_nom_k) - Hk*x_;
%         innovation(3:3:end) = mod(pi + innovation(3:3:end), 2*pi) - pi;
%         x = x_ + K*innovation;
%         P = (I-K*Hk)*P_;
%         inns(:,k) = innovation(1:min(end,3));
%     else
%         x = x_;
%         P = P_;
%         invS = inf;
%         innovation = inf;
%     end
%     
%     % store current values
%     xs(:,k) = x;
%     Ps(:,:,k) = P;
%     s(:,k) = 2*sqrt(diag(P));
%     invSs{k} = invS;
%     dy{k} = innovation;
% end
% invSs{1} = 0;
% dy{1} = 0;
% 
% % reconstruct full state estimate 
% xs = x_nom(t) + xs;
% end
% 
% 

% function [xs,Ps,s,invSs,ey,inns] = EKF(x0,P0,f,h,A,H,W,Q,R,t,ydata)
% opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
% ny = length(ydata);
% dt = t(2)-t(1);
% x = x0;
% P = P0;
% I = eye(length(x0));
% xs = zeros(4,ny); xs(:,1) = x;
% Ps = zeros(4,4,ny); Ps(:,:,1) = P;
% invSs = cell(1,ny);
% ey = cell(1,ny);
% s = zeros(4,ny); s(:,1) = 2*sqrt(diag(P));
% inns = zeros(3,ny);
% for k = 2:ny
%     % get F
%     Fk = I + dt*A(x);
%     
%     % prediction step
%     [~,x_] = ode45(f,[0 dt],x,opts);
%     x_ = x_(end,:)';
%     P_ = Fk*P*Fk' + W*Q*W';
% 
%     % measurement processing
%     if ~isempty(ydata{k})
%         i_in_view = ydata{k}(4,:);
%         yk = reshape(ydata{k}(1:3,:),[],1);
%     else
%         i_in_view = [];
%         yk = [];
%     end
%     Hs = cell(1,length(i_in_view)); Rs = Hs;
%     yk_est = zeros(3*length(i_in_view),1);
%     for j = 1:length(i_in_view)
%         Hs{j} = H(t(k),x,i_in_view(j));
%         Rs{j} = R;
%         yk_est((j-1)*3+(1:3)) = h(t(k),x_,i_in_view(j));
%     end
% 
%     % create augmented H and R matrices
%     Hk = vertcat(Hs{:});
%     Rk = blkdiag(Rs{:});
% 
%     % correction step
%     if ~isempty(i_in_view)
%         invS = inv(Hk*P_*Hk' + Rk);
%         K = P_*Hk'/(Hk*P_*Hk' + Rk);
%         innovation = yk - yk_est;
%         innovation(3:3:end) = mod(pi + innovation(3:3:end), 2*pi) - pi;
%         x = x_ + K*innovation;
%         P = (I-K*Hk)*P_;
%         inns(:,k) = innovation(1:min(end,3));
%         
%     else
%         x = x_;
%         P = P_;
%         invS = inf;
%         innovation = inf;
%     end
% 
%     % store current values
%     xs(:,k) = x;
%     Ps(:,:,k) = P;
%     s(:,k) = 2*sqrt(diag(P));
%     invSs{k} = invS;
%     ey{k} = innovation;
% end
% ey{1} = 0;
% invSs{1} = 0;
% end
% 
