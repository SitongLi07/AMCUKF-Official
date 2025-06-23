
% SOC Estimation from Ah with Zero Point Correction
% File: plot_soc_comparison.m

load('11-05-18_20.44 556_Charge2'); 
Q_nominal = 3.0; 

ah = meas.Ah;
if isrow(ah)
    ah = ah';
end

ah_zeroed = ah - ah(1);
Q_discharged_zeroed = ah_zeroed;
soc_corrected = 0.05 + Q_discharged_zeroed / Q_nominal;

figure;
plot(soc_corrected, '+-', 'DisplayName', 'Corrected SOC (Ah zeroed)');
xlabel('Sample Index');
ylabel('State of Charge (SOC)');
title('SOC Comparison: Zeroed-Ah Corrected');
grid on;
legend('show');
