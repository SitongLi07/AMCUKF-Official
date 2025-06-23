%Nonlinear Numerical System
close all;
clear all;
clc;
P0 = 5; % Initial noise covariance
Q = 1; % Process noise covariance
R = 9; % Measurement noise covariance
R1 = 100 * R;
tf = 100; % Final time
n_mc = 100;
h = inline('(x.^2)/20');
f = inline('0.5*x + 5*x./(1+x.^2) + 8*cos(0.4*t)','x','t');

randn('state',0);
x0 =sqrtm(P0)*randn(1); % Initial state value
sigmaF = 8;
nu = 5;

%% ---------- Select Noise Scenario ----------
% 1 = stationary‑t   noise
% 2 = nonstationary‑t noise
% 3 = Cauchy noise
sceneID = 1;                        
noiseFile = sprintf('dis%d.mat', sceneID);

%% ---------- Generate and Save Noise if File Not Found ----------
if exist(noiseFile, 'file')
    load(noiseFile, 'disturb_saved');           
    fprintf('Loaded fixed noise from  %s\n', noiseFile);
else
    fprintf('Generating noise for scene %d ...\n', sceneID);
    disturb_saved = zeros(1, tf);              

    switch sceneID
        case 1   % stationary t noise
            for t = 1:tf
                disturb_saved(t) = multivariate_t_noise(R1, nu, 1);
            end
        case 2   % non‑stationary t noise
            for t = 1:tf
                if   t <= 30
                    disturb_saved(t) = multivariate_t_noise(50*R, nu, 1);
                elseif t <= 70
                    disturb_saved(t) = multivariate_t_noise(100*R, nu, 1);
                else
                    disturb_saved(t) = multivariate_t_noise(30*R, nu, 1);
                end
            end
        case 3   % Cauchy noise
            disturb_saved = cauchy_noise(0, 5, 1, tf);
    end

    save(noiseFile, 'disturb_saved');           
    fprintf('Noise saved to  %s\n', noiseFile);
end


for j = 1 : n_mc
    j
    for t = 1:tf % Simulate the system
        if t == 1
            x(t) = feval(f,x0,t) + sqrtm(Q)*randn(1);
        else
            x(t) = feval(f,x(t-1),t) + 1*sqrtm(Q)*randn(1);% 
        end

%%-----------------------------------------------------------------------%%
%         disturb = (multivariate_t_noise(R1, nu, 1))'; % t noise 
%         y(t) = feval(h,x(t)) + disturb;
%%-----------------------------------------------------------------------%%
%         if t <= 30  % nostable noise 
%             disturb = (multivariate_t_noise(50*R, nu, 1))';
%         elseif t <= 70
%             disturb = (multivariate_t_noise(100*R, nu, 1))';
%         else 
%             disturb = (multivariate_t_noise(40*R, nu, 1))';
%         end
%         y(t) = feval(h,x(t)) + disturb;
%%-----------------------------------------------------------------------%%
%         disturb = cauchy_noise(0, 5, 1, 1); % cauchy_noise
%         y(t) = feval(h,x(t))+disturb;
%%-----------------------------------------------------------------------%%
        disturb        = disturb_saved(t);  
        y(t)           = feval(h, x(t)) + disturb;
        disturb_t(t) = disturb;
    end
    disturb_mc(j, :) = disturb_t;
    xTrue = x;
    xTrue_mc(j, :) = xTrue;
    % EKF
    tic
    xEKFhat = EKF(f, h, Q, R, x0, P0, y);
    elapsedTime = toc;
    xEKFrunningtime(1,j) = elapsedTime;
    xEKFhat = xEKFhat;
    xEKFhat_mc(j, :) = xEKFhat;
    xEKFrmse_mc(j, :) = xEKFhat - xTrue;

    tic
    xhat = UKF_A1(f, h, Q, R, x0, P0, y);
    elapsedTime = toc;
    xrunningtime(1,j) = elapsedTime;
    xhat = xhat;
    xhat_mc(j, :) = xhat;
    xrmse_mc(j, :) = xhat - xTrue;

    tic
    [xFhat, iterF_num,Lk_t] = FUKF_iter(f, h, Q, R, x0, P0, y,sigmaF);
    elapsedTime = toc;
    xFrunningtime(1,j) = elapsedTime;
    xFhat = xFhat;
    xFhat_mc(j, :) = xFhat;
    xFitem_mc(j,:) = iterF_num;
    xFrmse_mc(j, :) = xFhat - xTrue;
    Lk_mc(j,:) = Lk_t;

    tic
    [xEIhat, iterE_num,sigma_t]= EUKF_iter(f, h, Q, R, x0, P0, y);
    elapsedTime = toc;
    xErunningtime(1,j) = elapsedTime;
    xEIhat = xEIhat;
    xEIhat_mc(j, :) = xEIhat;
    xEitem_mc(j,:) = iterE_num;
    xEIrmse_mc(j, :) = xEIhat - xTrue;
    sigma_mc(j,:) = sigma_t;

    tic
    [xAhat, iterA_num,tau_t] = AUKF(f, h, Q, R, x0, P0, y);
    elapsedTime = toc;
    xArunningtime(1,j) = elapsedTime;
    xAhat = xAhat;
    xAhat_mc(j, :) = xAhat;
    xAitem_mc(j,:) = iterA_num;
    xArmse_mc(j, :) = xAhat - xTrue;
    tau_mc(j,:) = tau_t;
end
xEKFrmse_mean = mean(abs(xEKFrmse_mc), 1);
maxEKFrmse_error = mean(max(abs(xEKFrmse_mc), [], 2));
xrmse_mean = mean(abs(xrmse_mc), 1);
maxrmse_error = mean(max(abs(xrmse_mc), [], 2));
xFrmse_mean = mean(abs(xFrmse_mc), 1);
maxFrmse_error = mean(max(abs(xFrmse_mc), [], 2));
xFitem_mean = mean(xFitem_mc, 1);
xEIrmse_mean = mean(abs(xEIrmse_mc), 1);
maxEIrmse_error = mean(max(abs(xEIrmse_mc), [], 2));
xEitem_mean = mean(xEitem_mc, 1);
xArmse_mean = mean(abs(xArmse_mc), 1);
maxArmse_error = mean(max(abs(xArmse_mc), [], 2));
xAitem_mean = mean(xAitem_mc, 1);

% set(groot,'defaultTextInterpreter','latex')
% figure
% plot(1:tf, xTrue,'k', 1:tf,xEKFhat,'m--', 1:tf,xhat,'c:',1:tf,xFhat,'g-.*',1:tf,xEIhat,'bo:',1:tf,xAhat,'r-.', 'LineWidth', 1.5, 'MarkerSize', 3);
% xlabel('Time(s)');
% ylabel('State values(-)')
% legend('True states','EKF', 'UKF','FMCUKF with \sigma =8','EMCUKF','AMCUKF');
% % title('UKF Simulation by Yuan-Li Cai');
% % grid on;
% gcaall = gca;
% subAP1 = [40, -35, 20, 10];zoomAP1 = [22, -18, 5, 8];         
% zp = BaseZoom(gcaall, subAP1, zoomAP1);                                
% zp.run; 
set(groot,'defaultTextInterpreter','latex')
figure; hold on
plot(1:tf, xTrue, '-', 'Color', 'k', 'LineWidth', 1.5)
% % plot(1:tf, xEKFhat, '--', 'Color', [0, 0.447, 0.741], 'LineWidth', 1.5)
% % plot(1:tf, xhat, ':', 'Color', [0.85, 0.325, 0.098], 'LineWidth', 1.5)
% % plot(1:tf, xFhat, '-.*', 'Color', [0.929, 0.694, 0.125], 'LineWidth', 1.5, 'MarkerSize', 3)
% % plot(1:tf, xEIhat, 'o--', 'Color', [0.494, 0.184, 0.556], 'LineWidth', 1.5, 'MarkerSize', 3)
% % plot(1:tf, xAhat, '^', 'Color', [0.301, 0.745, 0.933], 'Linestyle', 'none', 'LineWidth', 1.5, 'MarkerSize', 3)
plot(1:tf, xEKFhat, '--', 'Color', [0, 0.45, 0.7], 'LineWidth', 1.5)           
plot(1:tf, xhat, '-.', 'Color', [0.9, 0.6, 0], 'LineWidth', 1.5)              
plot(1:tf, xEIhat, ':', 'Color', [0, 0.6, 0.5], 'LineWidth', 1.5)               
plot(1:tf, xFhat, '-s', 'Color', [0.8, 0.4, 0.6], 'LineWidth', 1.2, 'MarkerSize', 4, 'MarkerIndices', 1:6:tf, 'MarkerFaceColor', [0.8, 0.4, 0.6])
plot(1:tf, xAhat, '--^', 'Color', [0.6, 0.8, 0.2], 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerIndices', 1:5:tf, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.6, 0.8, 0.2]);



xlabel('Time(s)');
ylabel('State values(-)')
legend('True states','EKF', 'UKF','EMCUKF','FMCUKF(\sigma=8)','AMCUKF');
% axis([0, 100, -30, 60])
box on     
gcaall = gca;
% subAP1 = [45, -35, 15, 15];zoomAP1 = [21, -17, 6, 6];        
subAP1 = [33, -35, 15, 15];zoomAP1 = [21, -17, 6, 6];  
zp = BaseZoom(gcaall, subAP1, zoomAP1);                              
zp.run; 


set(groot,'defaultTextInterpreter','latex')
figure;hold on
% plot(1:tf, xEKFrmse_mean, '--', 'Color', [0, 0.447, 0.741], 'LineWidth', 1.5)
% plot(1:tf, xrmse_mean, ':', 'Color', [0.85, 0.325, 0.098], 'LineWidth', 1.5)
% plot(1:tf, xFrmse_mean, '-.*', 'Color', [0.929, 0.694, 0.125], 'LineWidth', 1.5, 'MarkerSize', 3)
% plot(1:tf, xEIrmse_mean, 'o--', 'Color', [0.494, 0.184, 0.556], 'LineWidth', 1.5, 'MarkerSize', 3)
% plot(1:tf, xArmse_mean, '^-', 'Color', [0.301, 0.745, 0.933], 'LineWidth', 1.5, 'MarkerSize', 3)
plot(1:tf, xEKFrmse_mean, '--', 'Color', [0, 0.45, 0.7], 'LineWidth', 1.5)          
plot(1:tf, xrmse_mean, '-.', 'Color', [0.9, 0.6, 0], 'LineWidth', 1.5)              
plot(1:tf, xEIrmse_mean, ':', 'Color', [0, 0.6, 0.5], 'LineWidth', 1.5)               
plot(1:tf, xFrmse_mean, '-s', 'Color', [0.8, 0.4, 0.6], 'LineWidth', 1.2, 'MarkerSize', 5, 'MarkerIndices', 1:6:tf, 'MarkerFaceColor', [0.8, 0.4, 0.6]) 
plot(1:tf, xArmse_mean, '--^', 'Color', [0.6, 0.8, 0.2], 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerIndices', 1:5:tf, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.6, 0.8, 0.2]);


xlabel('Time(s)');
ylabel('RMSE(-)')
legend('EKF', 'UKF','EMCUKF','FMCUKF(\sigma=8)','AMCUKF');
axis([0, 100, 0, 20]);
box on     
gcaall = gca;
subAP1 = [30, 14, 20, 4];zoomAP1 = [56, 0.4, 5, 0.9];
% subAP1 = [30, 14, 20, 5];zoomAP1 = [36, 0.75, 8, 1.25];
% subAP1 = [30, 14, 10, 4];zoomAP1 = [43.8, 0.75, 3, 1.2];        
zp = BaseZoom(gcaall, subAP1, zoomAP1);                                
zp.run;

figure
plot(1:tf, sqrt(tau_t), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3)
xlabel('Time (s)')
ylabel('Kernel Scale $\sigma_{k}$ (-)','Color','b')
hold on;

yyaxis right;
ylim ([-20, 20]);
for t = 1 : tf
  plot([t, t], [0, disturb_t(t)], 'k-')
  hold on;
end

hold off;


ylabel('Disturbance Noise(-)', 'Color', 'k')
legend('AMCUKF', 'Disturbance Noise')

% figure
% plot(1:tf,sqrt(sigma_nst),'m--');
% xlabel('Time');
% ylabel('sigma');

figure
subplot(3, 2, 1)
plot(1:tf, sigmaF*ones(1,tf), 'b-', 'LineWidth', 1.5, 'MarkerSize', 3)
xlabel(['Time (s)' char(10) '(a)'])
ylabel('$\sigma_{k}$')
legend('FMCUKF')
subplot(3, 2, 2)
plot(1:tf, Lk_t, 'k-.', 'LineWidth', 1.5, 'MarkerSize', 3)
xlabel(['Time (s)' char(10) '(b)'])
ylabel('$L_{k}$')
legend('FMCUKF')
axis([0, 100, 0, 1.4])
subplot(3, 2, 3)
plot(1:tf, sigma_t, 'b-', 'LineWidth', 1.5, 'MarkerSize', 3)
xlabel(['Time (s)' char(10) '(c)'])
ylabel('$\sigma_{k}$')
legend('EMCUKF')
subplot(3, 2, 4)
plot(1:tf, 0.6065*ones(1,tf), 'k-.', 'LineWidth', 1.5, 'MarkerSize', 3)
xlabel(['Time (s)' char(10) '(d)'])
ylabel('$L_{k}$')
legend('EMCUKF')
subplot(3, 2, 5)
plot(1:tf, sqrt(tau_t), 'b-', 'LineWidth', 1.5, 'MarkerSize', 3)
xlabel(['Time (s)' char(10) '(e)'])
ylabel('$\sigma_{k}$')
legend('AMCUKF')
subplot(3, 2, 6)
plot(1:tf, 1./tau_t, 'k-.', 'LineWidth', 1.5, 'MarkerSize', 3)
xlabel(['Time (s)' char(10) '(f)'])
ylabel('$L_{k}$')
legend('AMCUKF')
axis([0, 100, 0.002, 0.027])

fprintf(' EKF RMSE: %.4f MEE: %.4f AT: %.4f ms\n', mean(xEKFrmse_mean), maxEKFrmse_error, mean(xEKFrunningtime)*10);
fprintf(' UKF RMSE: %.4f MEE: %.4f AT: %.4f ms\n', mean(xrmse_mean), maxrmse_error, mean(xrunningtime)*10);
fprintf('FUKF RMSE: %.4f MEE: %.4f AT: %.4f ms ANI: %.2f\n', mean(xFrmse_mean), maxFrmse_error, mean(xFrunningtime)*10, mean(xFitem_mean));
fprintf('EUKF RMSE: %.4f MEE: %.4f AT: %.4f ms ANI: %.2f\n', mean(xEIrmse_mean), maxEIrmse_error, mean(xErunningtime)*10, mean(xEitem_mean));
fprintf('AUKF RMSE: %.4f MEE: %.4f AT: %.4f ms ANI: %.2f\n', mean(xArmse_mean), maxArmse_error, mean(xArunningtime)*10, mean(xAitem_mean));
