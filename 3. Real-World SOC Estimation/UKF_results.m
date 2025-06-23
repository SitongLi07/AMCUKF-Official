
% Real-World SOC Estimation

clear; clc;close all;

data = readtable('SOC_Ah_Voltage_Current.csv');
current = data.Current;
% voltage = load("voltage_noise.mat").voltage;

% voltage = data.Voltage;
% 
% tnoise = 0.16 * trnd(3, size(voltage, 1), 1);
% voltage = voltage + tnoise;
% 
% save('voltage_noise02.mat', 'voltage');

voltage = load("voltage_noise02.mat").voltage;

soc_true = data.SOC;
SOCnoise = sqrt(1e-7) * rand(size(soc_true));
soc_true = soc_true + SOCnoise;

soc0 = soc_true(1);  

X = [0    0.0039    6.7835    0.0010    0.0044    3.0230];

tic;
[SOC_est_UKF, V_est_UKF, P_trace, state_history] = simulate_ESC_UKF(X, current, voltage, soc0);
elapedTime = toc;
fprintf('UKF run time：%.4f s\n', elapedTime);

tic;
[SOC_est_EUKF, V_est_EUKF, P_trace, state_history,iter_num] = simulate_ESC_EUKF(X, current, voltage, soc0);
elapedTime = toc;
fprintf('EMCUKF run time：%.4f s\n', elapedTime);
fprintf('EMCUKF number of iterationsrun time：%.4f \n', mean(iter_num));

tic;
[SOC_est_FUKF, V_est_FUKF, P_trace, state_history,iter_num] = simulate_ESC_FUKF(X, current, voltage, soc0, 8);
elapedTime = toc;
fprintf('FMCUKF run time：%.4f s\n', elapedTime);
fprintf('FUKF number of iterations：%.4f \n', mean(iter_num));

tic;
[SOC_est_AUKF, V_est_AUKF, P_trace, state_history,iter_num] = simulate_ESC_AUKF(X, current, voltage, soc0);
elapedTime = toc;
fprintf('AMCUKF run time：%.4f s\n', elapedTime);
fprintf('AMCUKF number of iterations：%.4f \n', mean(iter_num));


error_vec_UKF = SOC_est_UKF - soc_true;
rmse = sqrt(mean(error_vec_UKF.^2));
maxUKFrmse_error = max(abs(error_vec_UKF));
fprintf('UKF SOC RMSE = %.6f\n', rmse);
fprintf('UKF SOC ME = %.6f\n', maxUKFrmse_error);


error_vec_EUKF = SOC_est_EUKF - soc_true;
rmse = sqrt(mean(error_vec_EUKF.^2));
maxEUKFrmse_error = max(abs(error_vec_EUKF));
fprintf('EMCUKF SOC RMSE = %.6f\n', rmse);
fprintf('EMCUKF SOC ME = %.6f\n', maxEUKFrmse_error);


error_vec_FUKF = SOC_est_FUKF - soc_true;
rmse = sqrt(mean(error_vec_FUKF.^2));
maxFUKFrmse_error = max(abs(error_vec_FUKF));
fprintf('FMCUKF SOC RMSE = %.6f\n', rmse);
fprintf('FMCUKF SOC ME = %.6f\n', maxFUKFrmse_error);


error_vec_AUKF = SOC_est_AUKF - soc_true;
rmse = sqrt(mean(error_vec_AUKF.^2));
maxAUKFrmse_error = max(abs(error_vec_AUKF));
fprintf('AMCUKF SOC RMSE = %.6f\n', rmse);
fprintf('AMCUKF SOC ME = %.6f\n', maxAUKFrmse_error);





set(groot,'defaultTextInterpreter','latex')
figure;
plot(soc_true, 'k', 'LineWidth', 1.5, 'DisplayName', 'True SOC'); hold on;
plot(SOC_est_UKF, '-.', 'Color', [0.9, 0.6, 0], 'LineWidth', 1.5,'DisplayName', 'UKF')               
plot(SOC_est_EUKF, ':', 'Color', [0, 0.6, 0.5], 'LineWidth', 1.5,'DisplayName', 'EMCUKF')               
plot(SOC_est_FUKF, '-s', 'Color', [0.8, 0.4, 0.6], 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerIndices', 1:6:length(SOC_est_FUKF), 'MarkerFaceColor', [0.8, 0.4, 0.6],'DisplayName', 'FMCUKF(\sigma=8)') 
plot(SOC_est_AUKF, '--^', 'Color', [0.6, 0.8, 0.2], 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerIndices', 1:5:length(SOC_est_AUKF), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.6, 0.8, 0.2],'DisplayName', 'AMCUKF');

xlabel('Times (min)'); 
ylabel('State of Charge, SOC (-)');
legend('show'); 
grid on;
axis([0, 88, 0, 1.3]);
gcaall = gca;
subAP1 = [10, 0.7, 15, 0.2];zoomAP1 = [21, 0.37, 3, 0.06];         
zp = BaseZoom(gcaall, subAP1, zoomAP1);                                
zp.run; 






rmse_seq_UKF = sqrt(error_vec_UKF.^2);
rmse_seq_EUKF = sqrt(error_vec_EUKF.^2);
rmse_seq_FUKF = sqrt(error_vec_FUKF.^2);
rmse_seq_AUKF = sqrt(error_vec_AUKF.^2);

set(groot,'defaultTextInterpreter','latex')
figure;hold on;

plot(rmse_seq_UKF, '-.', 'Color', [0.9, 0.6, 0], 'LineWidth', 1.5,'DisplayName', 'UKF')               
plot(rmse_seq_EUKF, ':', 'Color', [0, 0.6, 0.5], 'LineWidth', 1.5,'DisplayName', 'EMCUKF')               
plot(rmse_seq_FUKF, '-s', 'Color', [0.8, 0.4, 0.6], 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerIndices', 1:6:length(SOC_est_FUKF), 'MarkerFaceColor', [0.8, 0.4, 0.6],'DisplayName', 'FMCUKF(\sigma=8)') 
plot(rmse_seq_AUKF, '--^', 'Color', [0.6, 0.8, 0.2], 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerIndices', 1:5:length(SOC_est_AUKF), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.6, 0.8, 0.2],'DisplayName', 'AMCUKF');

xlabel('Times (min)');
ylabel('RMSE (-)');
legend('show');
grid on;
axis([0, 88, -0.005, 0.09]);
box on