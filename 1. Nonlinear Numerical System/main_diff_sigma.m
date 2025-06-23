% Generates kernel scale comparison plots (Fig. 3)
close all;
clear all;
clc;
P0 = 5; % Initial noise covariance
Q = 10; % Process noise covariance
R = 1; % Measurement noise covariance
tf = 100; % Final time
n_mc = 20;
sigma_bag = [0.1:0.3:1, 2:2:30];
h = inline('(x.^2)/20');
f = inline('0.5*x + 5*x./(1+x.^2) + 8*cos(0.4*t)','x','t');
randn('state',0);
x0 =sqrtm(P0)*randn(1); % Initial state value
% sigmaF = 1:n_mc;
% xTrue_mc = zeros();
for i = 1 : size(sigma_bag,2)
    sigmaF = sigma_bag(:,i);
    for j = 1 : n_mc

        for t = 1:tf % Simulate the system
            if t == 1
                x(t) = feval(f,x0,0) + sqrtm(Q)*randn(1);
            else
                x(t) = feval(f,x(t-1),t-1) + 1*sqrtm(Q)*randn(1);% + 0.5*sqrtm(2*Q)*randn(1);
            end
            disturb = 1*sqrtm(100*R)*randn(1); % large Gaussian noise
            % disturb = randi([-20,20]);
            % y(t) = feval(h,x(t)) + 1*sqrtm(R)*randn(1)+1*sqrtm(100*R)*randn(1);

            y(t) = feval(h,x(t)) + 1*sqrtm(R)*randn(1)+disturb;
            disturb_t(t) = disturb;
        end
        disturb_mc(j, :) = disturb_t;
        xTrue = [x0, x];
        xTrue_mc(j, :) = xTrue;

        % xhat = UKF_A1(f, h, Q, R, x0, P0, y);
        % xhat = [x0, xhat];
        % xhat_mc(j, :) = xhat;
        % xrmse_mc(j, :) = xhat - xTrue;

        % xEhat = EUKF(f, h, Q, R, x0, P0, y);
        % xEhat = [x0, xEhat];
        % xEhat_mc(j, :) = xEhat;
        % xErmse_mc(j, :) = xEhat - xTrue;
        [xFhat, iterF_num,Lk_t] = FUKF_iter(f, h, Q, R, x0, P0, y,sigmaF);
        xFhat = [x0, xFhat];
        xFhat_mc(j, :, i) = xFhat;
        xFitem_mc(j, :, i) = iterF_num;
        xFrmse_mc(j, :, i) = abs(xFhat - xTrue);
        Lk_mc(j,:, i) = Lk_t;

        % [xEIhat, iterE_num,sigma_t]= EUKF_iter(f, h, Q, R, x0, P0, y);
        % xEIhat = [x0, xEIhat];
        % xEIhat_mc(j, :) = xEIhat;
        % xEitem_mc(j,:) = iterE_num;
        % xEIrmse_mc(j, :) = xEIhat - xTrue;
        % sigma_mc(j,:) = sigma_t;
        %
        % [xAhat, iterA_num,tau_t] = AUKF(f, h, Q, R, x0, P0, y);
        % xAhat = [x0, xAhat];
        % xAhat_mc(j, :) = xAhat;
        % xAitem_mc(j,:) = iterA_num;
        % xArmse_mc(j, :) = xAhat - xTrue;
        % tau_mc(j,:) = tau_t;
    end
    xFrmse_mean_mc(i,:) =mean(xFrmse_mc(:,:, i),1);
end
xFrmse_mean_sigma = mean(xFrmse_mean_mc,2);
% xrmse_mean = mean(abs(xrmse_mc), 1);
% xErmse_mean = mean(abs(xErmse_mc), 1);
% xFrmse_mean = mean(abs(xFrmse_mc), 2);
% xFitem_mean = mean(xFitem_mc, 2);
% xEIrmse_mean = mean(abs(xEIrmse_mc), 1);
% xEitem_mean = mean(xEitem_mc, 1);
% xArmse_mean = mean(abs(xArmse_mc), 1);
% xAitem_mean = mean(xAitem_mc, 1);

set(groot,'defaultTextInterpreter','latex')
figure
plot(sigma_bag,xFrmse_mean_sigma','k-*')
xlabel('$\sigma_{k}$')
ylabel('ARMSE')
axis([0,30,4,9])