function [kfdx, P, phatm, dy, dy_pert, stdev, innov] = linearKalmanFilter(yStoreNoise, yStoreNominal, state0, P0, F, gamma, Q, R, H, deltaT)

kfdx = state0; P = P0; xhatm = []; phatm = []; dy = [];

% Ohm_k = dt*gamma; % Might Need to Check into This

stdev = 2*sqrt(diag(P0));

for i=1:(length(yStoreNoise)-1)
    xhatm(:,i+1) = F{i}*kfdx(i,:)';
    phatm(:,:,i+1) = F{i}*P(:,:,i)*F{i}' + Q;
    if(isempty(yStoreNoise{i+1})) || (isempty(yStoreNominal{i+1}))
        kfdx(i+1,:) = xhatm(:,i+1);
        P(:,:,i+1) = phatm(:,:,i+1);
        stdev(:,i+1) = 2*sqrt(diag(P(:,:,i+1)));
    else
        % yStoreNoise{i+1} = yStoreNoise{i+1}(1:3,:);
        if size(yStoreNominal{i+1},1) == 3
            currY = yStoreNoise{i+1}(1:3,1); currY = currY(:);
        else
            currY = yStoreNoise{i+1}(1:3,:); currY = currY(:);
        end
        [matSize, ~] = size(yStoreNoise{i+1});

        big_R = kron(eye(length(yStoreNominal{i+1})/3),R);
        K = phatm(:,:,i+1)*H{i+1}'*inv(H{i+1}*phatm(:,:,i+1)*H{i+1}' + big_R);

        if size(yStoreNominal{i+1},1) > 3 && size(currY,1) <= 3
            currY = [currY;0;0;0];
        end

        dy{i+1} = currY - yStoreNominal{i+1};
        dy_pert{i+1} = H{i+1}*xhatm(:,i+1);
        if matSize == 3 || matSize == 4
            % if matSize == 3
            dy{i+1}(3) = wrapToPi(dy{i+1}(3));
            dy_pert{i+1}(3) = wrapToPi(dy_pert{i+1}(3));
        else
            dy{i+1}(3) = wrapToPi(dy{i+1}(3));
            dy{i+1}(6) = wrapToPi(dy{i+1}(6));
            dy_pert{i+1}(3) = wrapToPi(dy_pert{i+1}(3));
            dy_pert{i+1}(6) = wrapToPi(dy_pert{i+1}(6));
        end

        innovation{i+1} = dy{i+1} - dy_pert{i+1};

        if matSize == 3 || matSize == 4
            innovation{i+1}(3) = wrapToPi(innovation{i+1}(3));
        else
            innovation{i+1}(6) = wrapToPi(innovation{i+1}(6));
        end

        kfdx(i+1,:) = (xhatm(:,i+1) + K*innovation{i+1});
        P(:,:,i+1) = (eye(4) - K*H{i+1})*phatm(:,:,i+1);
        innov{i+1} = innovation{i+1}(1:3);
        stdev(:,i+1) = 2*sqrt(diag(P(:,:,i+1)));
    end
end
end