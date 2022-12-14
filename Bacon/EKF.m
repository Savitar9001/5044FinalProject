function [xs,Ps,s,invSs,ey,inns] = EKF(x0,P0,f,h,A,H,Q,R,t,ydata)
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
ny = length(ydata);
dt = t(2)-t(1);
x_ = x0;
P_ = P0;
I = eye(length(x0));
xs = zeros(4,ny); 
Ps = zeros(4,4,ny);
invSs = cell(1,ny);
ey = cell(1,ny);
s = zeros(4,ny); 
inns = zeros(3,ny);
for k = 1:ny
    
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
        Hs{j} = H(t(k),x_,i_in_view(j));
        Rs{j} = R;
        yk_est((j-1)*3+(1:3)) = h(t(k),x_,i_in_view(j));
    end

    % create augmented H and R matrices
    Hk = vertcat(Hs{:});
    Rk = blkdiag(Rs{:});

    % correction step
    if ~isempty(i_in_view)
        invS = inv(Hk*P_*Hk' + Rk);
        K = P_*Hk'/(Hk*P_*Hk' + Rk);
        innovation = yk - yk_est;
        innovation(3:3:end) = mod(pi + innovation(3:3:end), 2*pi) - pi;
        x = x_ + K*innovation;
        P = (I-K*Hk)*P_;
        inns(:,k) = innovation(1:min(end,3));
        
    else
        x = x_;
        P = P_;
        invS = inf;
        innovation = NaN;
    end

    % get F
    Fk = I + dt*A(x);
    
    % prediction step
    [~,x_] = ode45(f,[0 dt],x,opts);
    x_ = x_(end,:)';
    P_ = Fk*P*Fk' + Q;


    % store current values
    xs(:,k) = x;
    Ps(:,:,k) = P;
    s(:,k) = 2*sqrt(diag(P));
    invSs{k} = invS;
    ey{k} = innovation;
end
end

