%UKF simulation example
close all;
clear all;
clc;

Q = 1*diag([1,1,1,1]); 
tf = 100; % Final time
n_mc = 20;

% h = inline('(x.^2)/20');
% f = inline('0.5*x + 5*x./(1+x.^2) + 8*cos(0.4*t)','x','t');
T = 0.5;
tf_num = tf / T;
sita = pi/60;
sita_T = sita * T;

F = [1, sin(sita_T)/sita, 0, -(1-cos(sita_T))/sita; 0, cos(sita_T), 0, -sin(sita_T); 0, (1-cos(sita_T))/sita, 1, sin(sita_T)/sita; 0, sin(sita_T), 0, cos(sita_T)];
G = [(T^2)/2, 0; T, 0; 0, (T^2)/2; 0, T];
sigma_vel = 5;
Q1 = diag([sigma_vel^2, sigma_vel^2]);
nu = 4;

Q = G * Q1 * G'; % Process noise covariance
R = diag([(pi/180)^2,10^2]); % Measurement noise covariance

randn('state',0);
x0 =[250;0;750;0];% + sqrtm(P0)*randn(4,1); % Initial state value
P0 = diag([50, 50, 50, 50]); % Initial noise covariance
sigmaF = 20;
n = size(x0, 1); m = size(R, 2);

sigma_bag = [0.1:0.5:20];
for i = 1 : size(sigma_bag,2)
    sigmaF = sigma_bag(:,i);
    for j = 1 : n_mc
        j
        for t = 1:tf_num % Simulate the system
            if t == 1
                % x(t) = feval(f,x0,0) + sqrtm(Q)*randn(1);
                x(:,t) = F * x0  + sigma_vel * G * randn(2, 1);
            else
                % x(t) = feval(f,x(t-1),t-1) + 1*sqrtm(Q)*randn(1);% + 0.5*sqrtm(2*Q)*randn(1);
                x(:, t) = F * x(:, t-1) +  sigma_vel * G * randn(2, 1);
            end
            % disturb = 1*sqrtm([9, 0; 0, 50]*R)*randn(m,1); % large Gaussian noise
            % disturb = randi([-1,20])*randn(2,1);
            disturb = (multivariate_t_noise(16*R, nu, 1))'; % t noise   
            % y(t) = feval(h,x(t)) + 1*sqrtm(R)*randn(1)+1*sqrtm(100*R)*randn(1);

            % y(t) = feval(h,x(t)) + 1*sqrtm(R)*randn(1)+disturb;
            y(:,t) = [atan2(x(1,t), x(3,t)); sqrt((x(1,t))^2+(x(3,t))^2)] + 1*sqrtm(R)*randn(m, 1)+disturb;
            disturb_t(:,t) = disturb;
        end
        disturb_mc(:, :, j) = disturb_t;
        xTrue = [x0, x];
        xTrue_mc(:, :, j) = xTrue;

        % xhat = UKF_A1(f, h, Q, R, x0, P0, y);
        % xhat = UKF_A1(F, Q, R, x0, P0, y);
        % xhat = [x0, xhat];
        % xhat_mc(:, :, j,i) = xhat;
        % xrmse_mc(:, :, j,i) = xhat - xTrue;
        % xrmseXY_mc(j, :,i) = sqrt((xrmse_mc(1,:, j,i)).^2 + (xrmse_mc(3,:, j,i)).^2);

        % xEhat = EUKF(f, h, Q, R, x0, P0, y);
        % xEhat = [x0, xEhat];
        % xEhat_mc(j, :) = xEhat;
        % xErmse_mc(j, :) = xEhat - xTrue;
        [xFhat, iterF_num,Lk_t] = FUKF_iter(F, Q, R, x0, P0, y,sigmaF);
        xFhat = [x0, xFhat];
        xFhat_mc(:, :, j, i) = xFhat;
        xFitem_mc(j, :, i) = iterF_num;
        xFrmse_mc(:, :, j, i) = xFhat - xTrue;
        xFrmseXY_mc(j, :, i) = sqrt((xFrmse_mc(1,:, j, i)).^2 + (xFrmse_mc(3,:, j, i)).^2);
        Lk_mc(j,:, i) = Lk_t;
        %
        % [xEIhat, iterE_num,sigma_t]= EUKF_iter(F, Q, R, x0, P0, y);
        % xEIhat = [x0, xEIhat];
        % xEIhat_mc(:, :, j) = xEIhat;
        % xEitem_mc(j,:) = iterE_num;
        % xEIrmse_mc(:, :, j) = xEIhat - xTrue;
        % xEIrmseXY_mc(j, :) = sqrt((xEIrmse_mc(1,:, j)).^2 + (xEIrmse_mc(3,:, j)).^2);
        % sigma_mc(j,:) = sigma_t;
        % %
        % [xAhat, iterA_num,tau_t] = AUKF(F, Q, R, x0, P0, y);
        % xAhat = [x0, xAhat];
        % xAhat_mc(:, :, j) = xAhat;
        % xAitem_mc(j,:) = iterA_num;
        % xArmse_mc(:, :, j) = xAhat - xTrue;
        % xArmseXY_mc(j, :) = sqrt((xArmse_mc(1,:, j)).^2 + (xArmse_mc(3,:, j)).^2);
        % tau_mc(j,:) = tau_t;
    end
    xFrmse_mean_mc(i,:) = mean(xFrmseXY_mc(:, :, i), 1);
end
xFrmse_mean = (mean(xFrmse_mean_mc, 2))';




set(groot,'defaultTextInterpreter','latex')

figure
plot(sigma_bag, xFrmse_mean, 'k-*')
xlabel('$\sigma_{k}$(-)')
ylabel('ARMSE(m)')
% axis([0,30,4,9])

% figure
% plot(xTrue(1,:), xTrue(3,:),'k', xhat(1,:),xhat(3,:),'m:', xFhat(1,:),xFhat(3,:),'g-.*', xEIhat(1,:),xEIhat(3,:),'bo:',...
%     xAhat(1,:),xAhat(3,:),'r--', 'LineWidth', 1.5, 'MarkerSize', 3);
% 
%     % 
% xlabel('Time');
% ylabel('State values')
% legend('True states','UKF', 'FMCUKF with \sigma =8','EMCUKF','AMCUKF');
% 
% figure
% plot(0:T:T*tf_num, xrmse_mean,'m:',0:T:T*tf_num, xFrmse_mean,'g-.*', 0:T:T*tf_num, xEIrmse_mean,'bo:',...
%     0:T:T*tf_num, xArmse_mean,'r--','LineWidth', 1.5, 'MarkerSize', 3)
% % 
% xlabel('Time');
% ylabel('RMSE')
% legend('UKF','FMCUKF with \sigma =8','EMCUKF','AMCUKF');
% 
% 
% % figure
% % plot(0:tf, xTrue,'k', 0:tf,xhat,'m:',0:tf,xFhat,'g-.*',0:tf,xEIhat,'bo:',0:tf,xAhat,'r--', 'LineWidth', 1.5, 'MarkerSize', 3);
% % xlabel('Time');
% % xlabel('Time');
% % ylabel('State values')
% % legend('True states','UKF','FMCUKF with \sigma =8','EMCUKF','AMCUKF');
% % % title('UKF Simulation by Yuan-Li Cai');
% % % grid on;
% % 
% % set(groot,'defaultTextInterpreter','latex')
% % figure
% % plot(0:tf, xrmse_mean,'m:',0:tf,xFrmse_mean,'g-.',0:tf, xEIrmse_mean,'bo:',0:tf, xArmse_mean,'r--','LineWidth', 1.5, 'MarkerSize', 3)
% % xlabel('Time');
% % ylabel('RMSE')
% % legend('UKF','FMCUKF with \sigma = 8','EMCUKF','AMCUKF');
% % 
% figure
% plot(1:1:T*tf_num, xFitem_mean(1, 2:2:end), 'g-.', 1:1:T*tf_num, xEitem_mean(1, 2:2:end), 'bo:', 1:1:T*tf_num, xAitem_mean(1, 2:2:end), 'r--','LineWidth', 1.5, 'MarkerSize', 3)
% xlabel('Time');
% ylabel('Number of iterations')
% legend('FMCUKF with \sigma = 8','EMCUKF','AMCUKF');
% % 
% % figure
% % plot(1:tf, sqrt(tau_t), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3)
% % xlabel('Time (s)')
% % ylabel('Kernel Scale $\sigma_{k}$','Color','b')
% % hold on;
% % 
% % yyaxis right;
% % ylim ([-20, 20]);
% % for t = 1 : tf
% %   plot([t, t], [0, disturb_t(t)], 'k-')
% %   hold on;
% % end
% % 
% % hold off;
% % 
% % ylabel('Disturbance Noise', 'Color', 'k')
% % legend('AMCUKF', 'Disturbance Noise')
% % 
% % figure
% % subplot(3, 2, 1)
% % plot(1:tf, sigmaF*ones(1,tf), 'r-', 'LineWidth', 1.5, 'MarkerSize', 3)
% % xlabel(['Time (s)' char(10) '(a)'])
% % ylabel('$\sigma_{k}$')
% % legend('FMCUKF')
% % subplot(3, 2, 2)
% % plot(1:tf, Lk_t, 'r-', 'LineWidth', 1.5, 'MarkerSize', 3)
% % xlabel(['Time (s)' char(10) '(b)'])
% % ylabel('$L_{k}$')
% % legend('FMCUKF')
% % axis([0, 100, 0, 1.4])
% % subplot(3, 2, 3)
% % plot(1:tf, sigma_t, 'r-', 'LineWidth', 1.5, 'MarkerSize', 3)
% % xlabel(['Time (s)' char(10) '(c)'])
% % ylabel('$\sigma_{k}$')
% % legend('EMCUKF')
% % subplot(3, 2, 4)
% % plot(1:tf, 0.6065*ones(1,tf), 'r-', 'LineWidth', 1.5, 'MarkerSize', 3)
% % xlabel(['Time (s)' char(10) '(d)'])
% % ylabel('$L_{k}$')
% % legend('EMCUKF')
% % subplot(3, 2, 5)
% % plot(1:tf, sqrt(tau_t), 'r-', 'LineWidth', 1.5, 'MarkerSize', 3)
% % xlabel(['Time (s)' char(10) '(e)'])
% % ylabel('$\sigma_{k}$')
% % legend('AMCUKF')
% % subplot(3, 2, 6)
% % plot(1:tf, 1./tau_t, 'r-', 'LineWidth', 1.5, 'MarkerSize', 3)
% % xlabel(['Time (s)' char(10) '(f)'])
% % ylabel('$L_{k}$')
% % legend('AMCUKF')
% % axis([0, 100, 0, 0.015])
% % 
% % xrms = mean(xrmse_mean)
% % xFrms = mean(xFrmse_mean)
% % % xErms = mean(xErmse_mean)
% % xEIrms = mean(xEIrmse_mean)
% % xArms = mean(xArmse_mean)
% % % rms = sqrt(sum(((xTrue - xhat).^2)/tf))