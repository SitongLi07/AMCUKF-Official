%UKF 
function [xhat] = UKF_A1(f, h, Q, R, x0, P0, y)
tf = size(y,2);
x =x0;
P = P0;
n = size(P, 2); m = size(R, 2);
for t = 1:tf
    % time update
    [sigm, wm, wc] = OnSigmaPoint(x, P); % generate sigma points
    sigmaX = 0.5*sigm + 5*sigm./(1+sigm.^2) + 8*cos(0.4*t);
    xhat_minus = sum(sigmaX.*repmat(wm,n,1));

    % P_minus = Q+sum(((sigmaX-xhat_minus)*(sigmaX -xhat_minus)').*repmat(wc,n,1));
    % sigmaY = feval(h, sigmaX);
    % y_minus = sum(sigmaY.*repmat(wm,n,1));
    % P_yy = R + sum(((sigmaY-y_minus)*(sigmaY-y_minus)').*repmat(wc,n,1));
    % P_xy = sum(((sigmaX-xhat_minus)*(sigmaY-y_minus)').*repmat(wc,n,1));

    P_minus = zeros(n, n);
    sigmaY = (sigmaX.^2)/20;
    y_minus = sum(sigmaY.*repmat(wm,n,1));
    P_yy = zeros(m,m);
    P_xy = zeros(n,m);
    for i = 1 : 2*n + 1
        P_minus = P_minus + wm(:, i) * (sigmaX(:, i)- xhat_minus) * (sigmaX(:, i)- xhat_minus)';
        P_yy = P_yy + wc(:, i) * (sigmaY(:, i)-y_minus) * (sigmaY(:, i)-y_minus)';
        P_xy = P_xy + wc(:, i) * (sigmaX(:, i)- xhat_minus) * (sigmaY(:, i)-y_minus)';
    end
    P_minus = P_minus + Q;
    P_yy = P_yy + R;
    

    Kg = P_xy*inv(P_yy);
    xhat(t) = xhat_minus + Kg*(y(t) - y_minus);
    x = xhat(t);
    P = P_minus - Kg*P_yy*Kg';
end