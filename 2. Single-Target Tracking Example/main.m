%Single-Target Tracking Example
close all;
clear all;
clc;

% Q = 1*diag([1,1,1,1]); 
tf = 100; % Final time
n_mc = 100;

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

Q = G * Q1 * G'; % Process noise covariance
R = diag([(pi/180)^2,5^2]); % Measurement noise covariance
nu = 4;

rng(0, 'twister'); 
x0 =[250;0;750;0];% + sqrtm(P0)*randn(4,1); % Initial state value
P0 = diag([50, 50, 50, 50]); % Initial noise covariance
sigmaF = 8;
n = size(x0, 1); m = size(R, 2);
% xTrue_mc = zeros();
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

%         disturb = 1*sqrtm([9, 0; 0, 50]*R)*randn(m,1); % large Gaussian noise
%         y(:,t) = [atan2(x(1,t), x(3,t)); sqrt((x(1,t))^2+(x(3,t))^2)] + 1*sqrtm(R)*randn(m, 1)+disturb;

        disturb = (multivariate_t_noise(16*R, nu, 1))'; % t noise        
        y(:,t) = [atan2(x(1,t), x(3,t)); sqrt((x(1,t))^2+(x(3,t))^2)] +disturb;

        disturb_t(:,t) = disturb;
    end
    disturb_mc(:, :, j) = disturb_t;
    xTrue = [x0, x];
    xTrue_mc(:, :, j) = xTrue;

    % EKF
    tic
    xEKFhat = EKF(F, Q, R, x0, P0, y);
    elapsedTime = toc;
    xEKFrunningtime(1,j) = elapsedTime;
    xEKFhat = [x0, xEKFhat];
    xEKFhat_mc(:, :, j) = xEKFhat;
    xEKFrmse_mc(:, :, j) = xEKFhat - xTrue;
    xEKFrmseXY_mc(j, :) = sqrt((xEKFrmse_mc(1,:, j)).^2 + (xEKFrmse_mc(3,:, j)).^2);

    tic
    xhat = UKF_A1(F, Q, R, x0, P0, y);
    elapsedTime = toc;
    xrunningtime(1,j) = elapsedTime;
    xhat = [x0, xhat];
    xhat_mc(:, :, j) = xhat;
    xrmse_mc(:, :, j) = xhat - xTrue;
    xrmseXY_mc(j, :) = sqrt((xrmse_mc(1,:, j)).^2 + (xrmse_mc(3,:, j)).^2);

    tic
    [xFhat, iterF_num,Lk_t] = FUKF_iter(F, Q, R, x0, P0, y,sigmaF);
    elapsedTime = toc;
    xFrunningtime(1,j) = elapsedTime;
    xFhat = [x0, xFhat];
    xFhat_mc(:, :, j) = xFhat;
    xFitem_mc(j, :) = iterF_num;
    xFrmse_mc(:, :, j) = xFhat - xTrue;
    xFrmseXY_mc(j, :) = sqrt((xFrmse_mc(1,:, j)).^2 + (xFrmse_mc(3,:, j)).^2);
    Lk_mc(j,:) = Lk_t;
    
    tic
    [xEIhat, iterE_num,sigma_t]= EUKF_iter(F, Q, R, x0, P0, y);
    elapsedTime = toc;
    xErunningtime(1,j) = elapsedTime;
    xEIhat = [x0, xEIhat];
    xEIhat_mc(:, :, j) = xEIhat;
    xEitem_mc(j,:) = iterE_num;
    xEIrmse_mc(:, :, j) = xEIhat - xTrue;
    xEIrmseXY_mc(j, :) = sqrt((xEIrmse_mc(1,:, j)).^2 + (xEIrmse_mc(3,:, j)).^2);
    sigma_mc(j,:) = sigma_t;
    
    tic
    [xAhat, iterA_num,tau_t] = AUKF(F, Q, R, x0, P0, y);
    elapsedTime = toc;
    xArunningtime(1,j) = elapsedTime;
    xAhat = [x0, xAhat];
    xAhat_mc(:, :, j) = xAhat;
    xAitem_mc(j,:) = iterA_num;
    xArmse_mc(:, :, j) = xAhat - xTrue;
    xArmseXY_mc(j, :) = sqrt((xArmse_mc(1,:, j)).^2 + (xArmse_mc(3,:, j)).^2);
    tau_mc(j,:) = tau_t;
end

xEKFrmse_mean = mean(xEKFrmseXY_mc, 1);
maxEKFrmse_error = mean(max(xEKFrmseXY_mc, [], 2));
fprintf(' EKF RMSE: %.4f MEE: %.4f AT: %.4f ms\n', mean(xEKFrmse_mean), maxEKFrmse_error, mean(xEKFrunningtime)*5);
xrmse_mean = mean(xrmseXY_mc, 1);
maxrmse_error = mean(max(xrmseXY_mc, [], 2));
fprintf(' UKF RMSE: %.4f MEE: %.4f AT: %.4f ms\n', mean(xrmse_mean), maxrmse_error, mean(xrunningtime)*5);
xFrmse_mean = mean(xFrmseXY_mc, 1);
xFitem_mean = mean(xFitem_mc, 1);
maxFrmse_error = mean(max(xFrmseXY_mc, [], 2));
fprintf('FUKF RMSE: %.4f MEE: %.4f AT: %.4f ms ANI: %.2f\n', mean(xFrmse_mean), maxFrmse_error, mean(xFrunningtime)*5, mean(xFitem_mean));
xEIrmse_mean = mean(xEIrmseXY_mc, 1);
xEitem_mean = mean(xEitem_mc, 1);
maxEIrmse_error = mean(max(xEIrmseXY_mc, [], 2));
fprintf('EUKF RMSE: %.4f MEE: %.4f AT: %.4f ms ANI: %.2f\n', mean(xEIrmse_mean), maxEIrmse_error, mean(xErunningtime)*5, mean(xEitem_mean));
xArmse_mean = mean(xArmseXY_mc, 1);
xAitem_mean = mean(xAitem_mc, 1);
maxArmse_error = mean(max(xArmseXY_mc, [], 2));
fprintf('AUKF RMSE: %.4f MEE: %.4f AT: %.4f ms ANI: %.2f\n', mean(xArmse_mean), maxArmse_error, mean(xArunningtime)*5, mean(xAitem_mean));

set(groot,'defaultTextInterpreter','latex')

figure
% plot(xTrue(1,:), xTrue(3,:),'k', xEKFhat(1,:),xEKFhat(3,:),'k^', xhat(1,:),xhat(3,:),'m:', xFhat(1,:),xFhat(3,:),'g-.*', xEIhat(1,:),xEIhat(3,:),'bo:',...
%     xAhat(1,:),xAhat(3,:),'r--', 'LineWidth', 1.5, 'MarkerSize', 3);
plot(xTrue(1,:), xTrue(3,:), '-', 'Color', 'k', 'LineWidth', 1.5); hold on;

plot(xEKFhat(1,:), xEKFhat(3,:), '--', 'Color', [0, 0.45, 0.7], 'LineWidth', 1.8); 
plot(xhat(1,:),   xhat(3,:),   '-.', 'Color', [0.9, 0.6, 0],   'LineWidth', 1.8); 
plot(xEIhat(1,:), xEIhat(3,:), ':',  'Color', [0, 0.6, 0.5],   'LineWidth', 1.8); 
plot(xFhat(1,:),  xFhat(3,:),  '-s', 'Color', [0.8, 0.4, 0.6], 'LineWidth', 1.4, ...
     'MarkerSize', 3, 'MarkerIndices', 1:6:length(xFhat), 'MarkerFaceColor', [0.8, 0.4, 0.6]); 
plot(xAhat(1,:),  xAhat(3,:),  '--^', 'Color', [0.6, 0.8, 0.2], 'LineWidth', 1.5, ...
     'MarkerSize', 3, 'MarkerIndices', 1:5:length(xAhat), 'MarkerFaceColor', 'none', ...
     'MarkerEdgeColor', [0.6, 0.8, 0.2]);  

    
xlabel('X (m)');
ylabel('Y (m)')
legend('True states','EKF', 'UKF','EMCUKF', 'FMCUKF (\sigma =8)','AMCUKF');
gcaall = gca;
% subAP1 = [45, -35, 15, 15];zoomAP1 = [21, -17, 6, 6];         
subAP1 = [500, -300, 750, 700];zoomAP1 = [625, -1475, 175, 125];  
zp = BaseZoom(gcaall, subAP1, zoomAP1);                              
zp.run; 


figure
% plot(0:T:T*tf_num, xEKFrmse_mean,'k^',0:T:T*tf_num, xrmse_mean,'m:',0:T:T*tf_num, xFrmse_mean,'g-.*', 0:T:T*tf_num, xEIrmse_mean,'bo:',...
%     0:T:T*tf_num, xArmse_mean,'r--','LineWidth', 1.5, 'MarkerSize', 3)
plot(0:T:T*tf_num, xEKFrmse_mean, '--', 'Color', [0, 0.45, 0.7], 'LineWidth', 1.5); hold on;   
plot(0:T:T*tf_num, xrmse_mean, '-.', 'Color', [0.9, 0.6, 0], 'LineWidth', 1.5);                
plot(0:T:T*tf_num, xEIrmse_mean, ':', 'Color', [0, 0.6, 0.5], 'LineWidth', 1.5);               

plot(0:T:T*tf_num, xFrmse_mean, '-s', 'Color', [0.8, 0.4, 0.6], 'LineWidth', 1.2, ...
    'MarkerSize', 4, 'MarkerIndices', 1:6:length(xFrmse_mean), 'MarkerFaceColor', [0.8, 0.4, 0.6]); 

plot(0:T:T*tf_num, xArmse_mean, '--^', 'Color', [0.6, 0.8, 0.2], 'LineWidth', 1.5, ...
    'MarkerSize', 5, 'MarkerIndices', 1:5:length(xArmse_mean), 'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', [0.6, 0.8, 0.2]); 

xlabel('Time(s)');
ylabel('RMSE(m)')
legend('EKF', 'UKF','EMCUKF', 'FMCUKF (\sigma =8)','AMCUKF');

% figure
% plot(0:tf, xTrue,'k', 0:tf,xhat,'m:',0:tf,xFhat,'g-.*',0:tf,xEIhat,'bo:',0:tf,xAhat,'r--', 'LineWidth', 1.5, 'MarkerSize', 3);
% xlabel('Time');
% xlabel('Time');
% ylabel('State values')
% legend('True states','UKF','FMCUKF with \sigma =8','EMCUKF','AMCUKF');
% % title('UKF Simulation by Yuan-Li Cai');
% % grid on;
% 
% set(groot,'defaultTextInterpreter','latex')
% figure
% plot(0:tf, xrmse_mean,'m:',0:tf,xFrmse_mean,'g-.',0:tf, xEIrmse_mean,'bo:',0:tf, xArmse_mean,'r--','LineWidth', 1.5, 'MarkerSize', 3)
% xlabel('Time');
% ylabel('RMSE')
% legend('UKF','FMCUKF with \sigma = 8','EMCUKF','AMCUKF');
% 
figure
plot(1:1:T*tf_num, xFitem_mean(1, 2:2:end), 'g-.', 1:1:T*tf_num, xEitem_mean(1, 2:2:end), 'bo:', 1:1:T*tf_num, xAitem_mean(1, 2:2:end), 'r--','LineWidth', 1.5, 'MarkerSize', 3)
xlabel('Time');
ylabel('Number of iterations')
legend('FMCUKF with \sigma = 8','EMCUKF','AMCUKF');
% 
% figure
% plot(1:tf, sqrt(tau_t), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3)
% xlabel('Time (s)')
% ylabel('Kernel Scale $\sigma_{k}$','Color','b')
% hold on;
% 
% yyaxis right;
% ylim ([-20, 20]);
% for t = 1 : tf
%   plot([t, t], [0, disturb_t(t)], 'k-')
%   hold on;
% end
% 
% hold off;
% 
% ylabel('Disturbance Noise', 'Color', 'k')
% legend('AMCUKF', 'Disturbance Noise')
% 
% figure
% subplot(3, 2, 1)
% plot(1:tf, sigmaF*ones(1,tf), 'r-', 'LineWidth', 1.5, 'MarkerSize', 3)
% xlabel(['Time (s)' char(10) '(a)'])
% ylabel('$\sigma_{k}$')
% legend('FMCUKF')
% subplot(3, 2, 2)
% plot(1:tf, Lk_t, 'r-', 'LineWidth', 1.5, 'MarkerSize', 3)
% xlabel(['Time (s)' char(10) '(b)'])
% ylabel('$L_{k}$')
% legend('FMCUKF')
% axis([0, 100, 0, 1.4])
% subplot(3, 2, 3)
% plot(1:tf, sigma_t, 'r-', 'LineWidth', 1.5, 'MarkerSize', 3)
% xlabel(['Time (s)' char(10) '(c)'])
% ylabel('$\sigma_{k}$')
% legend('EMCUKF')
% subplot(3, 2, 4)
% plot(1:tf, 0.6065*ones(1,tf), 'r-', 'LineWidth', 1.5, 'MarkerSize', 3)
% xlabel(['Time (s)' char(10) '(d)'])
% ylabel('$L_{k}$')
% legend('EMCUKF')
% subplot(3, 2, 5)
% plot(1:tf, sqrt(tau_t), 'r-', 'LineWidth', 1.5, 'MarkerSize', 3)
% xlabel(['Time (s)' char(10) '(e)'])
% ylabel('$\sigma_{k}$')
% legend('AMCUKF')
% subplot(3, 2, 6)
% plot(1:tf, 1./tau_t, 'r-', 'LineWidth', 1.5, 'MarkerSize', 3)
% xlabel(['Time (s)' char(10) '(f)'])
% ylabel('$L_{k}$')
% legend('AMCUKF')
% axis([0, 100, 0, 0.015])
% 
% xrms = mean(xrmse_mean)
% xFrms = mean(xFrmse_mean)
% % xErms = mean(xErmse_mean)
% xEIrms = mean(xEIrmse_mean)
% xArms = mean(xArmse_mean)
% % rms = sqrt(sum(((xTrue - xhat).^2)/tf))