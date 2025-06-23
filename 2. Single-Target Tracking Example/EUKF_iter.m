%EMCUKF
function [xhat, iter_num,sigma_t] = EUKF_iter(F, Q, R, x0, P0, y)
tf = size(y,2);
x =x0;
P = P0;
n = size(P, 2); m = size(R, 2);
% sigma = 5;
for t = 1:tf
    % time update
    [sigm, wm, wc] = OnSigmaPoint(x, P); % generate sigma points
    for i = 1 : size(sigm, 2)
        sigmaX(:,i) = F*sigm(:, i);
    end
    xhat_minus = sum(sigmaX.*repmat(wm,n,1),2);
    
    % P_minus = Q+sum(((sigmaX-xhat_minus)*(sigmaX -xhat_minus)').*repmat(wc,n,1));
    % sigmaY = feval(h, sigmaX);
    % y_minus = sum(sigmaY.*repmat(wm,n,1));
    % P_yy = R + sum(((sigmaY-y_minus)*(sigmaY-y_minus)').*repmat(wc,n,1));
    % P_xy = sum(((sigmaX-xhat_minus)*(sigmaY-y_minus)').*repmat(wc,n,1));
 
    P_minus = zeros(n, n);
    % sigmaY = feval(h, sigmaX);
    % y_minus = sum(sigmaY.*repmat(wm,n,1));
    for i = 1 : size(sigmaX,2)
        sigmaY(:, i) = [atan2(sigmaX(1,i), sigmaX(3,i)); sqrt((sigmaX(1,i))^2+(sigmaX(3,i))^2)];
    end
    y_minus = sum(sigmaY.*repmat(wm,m,1), 2);
    % P_yy = zeros(m,m);
    P_xy = zeros(n,m);
    for i = 1 : 2*n + 1
        P_minus = P_minus + wc(:, i) * (sigmaX(:, i)- xhat_minus) * (sigmaX(:, i)- xhat_minus)';
        % P_yy = P_yy + wc(:, i) * (sigmaY(:, i)-y_minus) * (sigmaY(:, i)-y_minus)';
        P_xy = P_xy + wc(:, i) * (sigmaX(:, i)- xhat_minus) * (sigmaY(:, i)-y_minus)';
    end
    P_minus = P_minus + Q;
    % P_yy = P_yy + R;
    
    
    x = xhat_minus;
    iter = 0;
    while(iter < 5) 
        xtemp = x;
        % y_minus = feval(h, x);
        y_err = y(:, t) - y_minus;
        sigma = y_err' * inv(R) * y_err;
        % Lk = exp(-0.5 / sigma);
        Lk = exp(-0.5);
       
        % P_yy_L = R + Lk*sum(((sigmaY-y_minus)*(sigmaY-y_minus)').*repmat(wc,n,1));
       
        P_yy = zeros(m,m);
        for i = 1 : 2*n + 1
            P_yy = P_yy + wc(:, i) * (sigmaY(:, i)-y_minus) * (sigmaY(:, i)-y_minus)';
        end
        P_yy_L = R + Lk*P_yy;
        
        Kg = Lk*P_xy*inv(P_yy_L);
        x = xhat_minus + Kg*(y(:,t) - y_minus);
        P = P_minus - Kg*P_yy*Kg';
        iter = iter + 1;
        if (abs(x - xtemp) / abs(xtemp) <= 1e-3)
            break;
        end
    end
    xhat(:,t) = x;
    iter_num(t) = iter;
    sigma_t(t) = sigma;
end