
% File: identify_ESC_parameters.m

clear; clc;

data = readtable('SOC_Ah_Voltage_Current.csv');
current = data.Current;
voltage = data.Voltage;
soc0 = data.SOC(1); 

% X = [M0, M, gamma, R0, R1, Cn]
lb = [0, 0, 0, 0.001, 0.001, 2.0];     
ub = [0.2, 0.2, 10, 0.1, 0.1, 4.0];    

x0 = [0.01, 0.01, 1.0, 0.01, 0.01, 3.0];


opts = optimoptions('particleswarm', 'Display', 'iter', 'SwarmSize', 30, 'MaxIterations', 100);

obj_fun = @(X) esc_objective_function(X, current, voltage, soc0);
[optX, fval] = particleswarm(obj_fun, 6, lb, ub, opts);


disp('Optimal parameters X = [M0, M, gamma, R0, R1, Cn]:');
disp(optX);
disp(['Minimum error J = ', num2str(fval)]);

[V_est, ~] = simulate_ESC_model(optX, current, soc0);
figure;
plot(voltage, 'k', 'DisplayName', 'Measured Voltage'); hold on;
plot(V_est, 'r--', 'DisplayName', 'ESC Model Estimate');
xlabel('Sample Index');
ylabel('Voltage (V)');
title('ESC Model Fitting Result');
legend('show');
grid on;
