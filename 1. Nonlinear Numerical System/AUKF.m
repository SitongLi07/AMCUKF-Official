%AMCUKF 
function [xhat,iter_num,tau_t] = AUKF(f, h, Q, R, x0, P0, y)
tf = size(y,2);
aalpha0 = 3; bbeta0 = 3; rhoA = 0.9; rhoB = 0.9;
aalpha = aalpha0; bbeta = bbeta0;
x =x0;
P = P0;
n = size(P, 2); m = size(R, 2);

for t = 1:tf
    % time update
    [sigm, wm, wc] = OnSigmaPoint(x, P); % generate sigma points
    sigmaX = 0.5*sigm + 5*sigm./(1+sigm.^2) + 8*cos(0.4*t);
    xhat_minus = sum(sigmaX.*repmat(wm,n,1));
    aalpha_minus = rhoA * aalpha; bbeta_minus = rhoB * bbeta;
    
    % P_minus = Q+sum(((sigmaX-xhat_minus)*(sigmaX -xhat_minus)').*repmat(wc,n,1));
    % sigmaY = feval(h, sigmaX);
    % y_minus = sum(sigmaY.*repmat(wm,n,1));
    % P_yy = R + sum(((sigmaY-y_minus)*(sigmaY-y_minus)').*repmat(wc,n,1));
    % P_xy = sum(((sigmaX-xhat_minus)*(sigmaY-y_minus)').*repmat(wc,n,1));
   
    P_minus = zeros(n, n);
    sigmaY = (sigmaX.^2)/20;
    y_minus = sum(sigmaY.*repmat(wm,n,1));
    % P_yy = zeros(m,m);
    P_xy = zeros(n,m);
    for i = 1 : 2*n + 1
        P_minus = P_minus + wm(:, i) * (sigmaX(:, i)- xhat_minus) * (sigmaX(:, i)- xhat_minus)';
        % P_yy = P_yy + wc(:, i) * (sigmaY(:, i)-y_minus) * (sigmaY(:, i)-y_minus)';
        P_xy = P_xy + wc(:, i) * (sigmaX(:, i)- xhat_minus) * (sigmaY(:, i)-y_minus)';
    end
    P_minus = P_minus + Q;
   
    
    x = xhat_minus; P = P_minus; aalpha =m/2 + aalpha_minus; bbeta = bbeta_minus;
    iter = 0;
    while(iter < 5)
        % temp = bbeta;
        xtemp = x;
        tau = bbeta / (aalpha - 1);
        % y_minus = feval(h, x);
      
        % P_yy_tau = tau*R + sum(((sigmaY-y_minus)*(sigmaY-y_minus)').*repmat(wc,n,1));
        % P_xy_tau = sum(((sigmaX-xhat_minus)*(sigmaY-y_minus)').*repmat(wc,n,1));

     
        P_yy = zeros(m,m);
        for i = 1 : 2*n + 1
            P_yy = P_yy + wc(:, i) * (sigmaY(:, i)-y_minus) * (sigmaY(:, i)-y_minus)';
        end
        P_yy_tau = tau*R + P_yy;
      
        Kg = P_xy*inv(P_yy_tau);
        x = xhat_minus + Kg*(y(t) - y_minus);
        P = P_minus - P_xy*inv(P_yy_tau)*P_xy';
        y_err = y(t) - (x.^2)/20;
        bbeta = bbeta_minus + 0.5*(y_err)' * inv(R) * (y_err);
        iter = iter + 1;
        % if (abs(temp - bbeta)/abs(temp) <= 1e-2)
        %     break;
        % end
        if (abs(x - xtemp) / abs(xtemp) <= 1e-3)
            break;
        end
    end
    xhat(t) = x;
    iter_num(t) = iter;
    tau_t(t) = bbeta / (aalpha -1);
end