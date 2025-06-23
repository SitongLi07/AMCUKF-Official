
% OCV-SOC 6th Order Polynomial Fitting
% File: fit_ocv_soc_poly6.m

data = readtable('SOC_Ah_Voltage_Current.csv');
soc = data.SOC;
voltage = data.Voltage;

coeffs = polyfit(soc, voltage, 6);

disp('6th-order polynomial coefficients (highest to lowest order):');
disp(coeffs');

soc_lin = linspace(0, 1, 200);
voltage_fit = polyval(coeffs, soc_lin);

figure;
plot(soc, voltage, '+', 'DisplayName', 'Raw Data', 'MarkerSize', 4); hold on;
plot(soc_lin, voltage_fit, '-', 'LineWidth', 2, 'DisplayName', '6th-Order Polynomial Fit');
xlabel('SOC');
ylabel('Voltage (V)');
title('OCV vs SOC - 6th Order Polynomial Fit');
grid on;
legend('show');
