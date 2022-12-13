function [kfdx, P, phatm, dy, dy_pert, stdev, innov] = linearKalmanFilter(y_store_noise, y_store_nom, dx, P_0, F, gamma, Q, R, H, dt)
    
    kfdx = dx; P = P_0; xhatm = []; phatm = []; dy = [];
    
    % Ohm_k = dt*gamma; % might need to use this
    
    stdev = 2*sqrt(diag(P_0));
    
    for i=1:(length(y_store_noise)-1)
       [mat_size, ~] = size(y_store_noise{i+1});
       
       xhatm(:,i+1) = F{i}*kfdx(i,:)';
       phatm(:,:,i+1) = F{i}*P(:,:,i)*F{i}' + Q;
       if(isempty(y_store_noise{i+1}))
          kfdx(i+1,:) = xhatm(:,i+1);
          P(:,:,i+1) = phatm(:,:,i+1);
          stdev(:,i+1) = 2*sqrt(diag(P(:,:,i+1)));
       else
           big_R = kron(eye(length(y_store_nom{i+1})/3),R);
           K = phatm(:,:,i+1)*H{i+1}'*inv(H{i+1}*phatm(:,:,i+1)*H{i+1}' + big_R);
           
           dy{i+1} = y_store_noise{i+1} - y_store_nom{i+1};
           dy_pert{i+1} = H{i+1}*xhatm(:,i+1);
           if mat_size == 3
               dy{i+1}(3) = wrapToPi(dy{i+1}(3));
               dy_pert{i+1}(3) = wrapToPi(dy_pert{i+1}(3));
           else
               dy{i+1}(3) = wrapToPi(dy{i+1}(3));
               dy{i+1}(6) = wrapToPi(dy{i+1}(6));
               dy_pert{i+1}(3) = wrapToPi(dy_pert{i+1}(3));
               dy_pert{i+1}(6) = wrapToPi(dy_pert{i+1}(6));
           end
        
      innovation{i+1} = dy{i+1} - dy_pert{i+1};
      
      if mat_size == 3
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