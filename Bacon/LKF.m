function [xs,Ps,s,invSs,dy,inns] = LKF(x0,P0,x_nom,h,A,H,W,Q,R,t,ydata)
% initial conditions
dt = t(2)-t(1);
ny = length(ydata);
x_ = x0;
P_ = P0;
xs = zeros(4,ny); 
Ps = zeros(4,4,ny);
invSs = cell(1,ny);
dy = cell(1,ny);
s = zeros(4,ny); 
inns = zeros(3,ny);
I = eye(length(x0));


for k = 1:ny
    % get model parameters
    x_nom_k = x_nom(t(k));
    
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
    
    
    % correct
    if ~isempty(i_in_view)
        invS = inv(Hk*P_*Hk' + Rk);
        K = P_*Hk'/(Hk*P_*Hk' + Rk);
        innovation = (yk-y_nom_k) - Hk*x_;
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

    Fk = I + dt*A(x_nom_k);
    % predict
    x_ = Fk*x;
    P_ = Fk*P*Fk' + W*Q*W';
    
    % store current values
    xs(:,k) = x;
    Ps(:,:,k) = P;
    s(:,k) = 2*sqrt(diag(P));
    invSs{k} = invS;
    dy{k} = innovation;
end


% reconstruct full state estimate 
xs = x_nom(t) + xs;
end
