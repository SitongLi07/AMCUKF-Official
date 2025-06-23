%UKF algorithm 1
function [xhat] = UKF_A1(F, Q, R, x0, P0, y)
tf = size(y,2);
x =x0;
P = P0;
n = size(P, 2); m = size(R, 2);
for t = 1:tf
    % time update
    [sigm, wm, wc] = OnSigmaPoint(x, P); % generate sigma points
    % sigmaX = feval(f, sigm, t);
    
    for i = 1 : size(sigm, 2)
        sigmaX(:,i) = F*sigm(:, i);
    end
    xhat_minus = sum(sigmaX.*repmat(wm,n,1),2);
    
%% for 循环代替 cai方法
    % P_minus = Q+sum(((sigmaX-xhat_minus)*(sigmaX -xhat_minus)').*repmat(wc,n,1));
    % sigmaY = feval(h, sigmaX);
    % y_minus = sum(sigmaY.*repmat(wm,n,1));
    % P_yy = R + sum(((sigmaY-y_minus)*(sigmaY-y_minus)').*repmat(wc,n,1));
    % P_xy = sum(((sigmaX-xhat_minus)*(sigmaY-y_minus)').*repmat(wc,n,1));
%%  更正方法
    P_minus = zeros(n, n);
    % sigmaY = feval(h, sigmaX);
    for i = 1 : size(sigmaX,2)
        sigmaY(:, i) = [atan2(sigmaX(1,i), sigmaX(3,i)); sqrt((sigmaX(1,i))^2+(sigmaX(3,i))^2)];
    end
    y_minus = sum(sigmaY.*repmat(wm,m,1), 2);
    P_yy = zeros(m,m);
    P_xy = zeros(n,m);
    for i = 1 : 2*n + 1
        P_minus = P_minus + wc(:, i) * (sigmaX(:, i)- xhat_minus) * (sigmaX(:, i)- xhat_minus)';
        P_yy = P_yy + wc(:, i) * (sigmaY(:, i)-y_minus) * (sigmaY(:, i)-y_minus)';
        P_xy = P_xy + wc(:, i) * (sigmaX(:, i)- xhat_minus) * (sigmaY(:, i)-y_minus)';
    end
    P_minus = P_minus + Q;
    P_yy = P_yy + R;
    %%

    Kg = P_xy*inv(P_yy);
    xhat(:,t) = xhat_minus + Kg*(y(:,t) - y_minus);
    x = xhat(:,t);
    P = P_minus - Kg*P_yy*Kg';
end