function xhat = EKF(F, H, Q, R, x0, P0, y)

tf = size(y,2);
x =x0;
P = P0;
n = size(P, 2); m = size(R, 2);

for t = 1:tf
    xpre = 0.5*x + 5*x./(1+x.^2) + 8*cos(0.4*t);
    F_jacx = 0.5 + 5*(1 - x.^2)./(1 + x.^2).^2;
    
    Ppre = F_jacx * P * F_jacx' + Q;
    H_jacxp = xpre / 10;

    zPre = (xpre.^2)/20;

    Pzz = H_jacxp * Ppre * H_jacxp' + R;
    Pxz = Ppre * H_jacxp';
    K = Pxz * inv(Pzz);

    x = xpre + K * (y(:, t) - zPre);
    P = Ppre - K * Pzz * K';

    xhat(:, t) = x;

end