function xhat = EKF(F, Q, R, x0, P0, y)

tf = size(y,2);
x =x0;
P = P0;
n = size(P, 2); m = size(R, 2);

for t = 1:tf
    xpre = F * x;
    Ppre = F * P * F' + Q;

    zPre = [atan2(xpre(1, :), xpre(3, :)); sqrt( (xpre(1, :))^2 + (xpre(3, :))^2 )];

    mag = (xpre(1, :))^2 + (xpre(3, :))^2;
    sqrt_mag = sqrt(mag);

    H = [xpre(3, :) / mag, 0, -xpre(1, :) / mag, 0;
         xpre(1, :) / sqrt_mag, 0, xpre(3, :) / sqrt_mag, 0];

    Pzz = H * Ppre * H' + R;
    Pxz = Ppre * H';
    K = Pxz * inv(Pzz);

    x = xpre + K * (y(:, t) - zPre);
    P = Ppre - K * Pzz * K';

    xhat(:, t) = x;



end