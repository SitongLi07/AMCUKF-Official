% Export SOC, Ah, Voltage, and Current to CSV 

load('11-05-18_20.44 556_Charge2.mat');
Q_nominal = 3.0;

Ah = meas.Ah;
Voltage = meas.Voltage;
Current = meas.Current;

if isrow(Ah), Ah = Ah'; end
if isrow(Voltage), Voltage = Voltage'; end
if isrow(Current), Current = Current'; end

idx_charge = Current > 0;
Ah = Ah(idx_charge);
Voltage = Voltage(idx_charge);
Current = Current(idx_charge);

% tnoise = 0.2 * trnd(3, size(Voltage, 1), 1);
% Voltage = Voltage + tnoise;

ah_zeroed = Ah - Ah(1);
Q_discharged_zeroed = ah_zeroed;
soc_corrected = 0.05 + Q_discharged_zeroed / Q_nominal;

T = table(soc_corrected, ah_zeroed, Voltage, Current, ...
    'VariableNames', {'SOC', 'Ah', 'Voltage', 'Current'});

writetable(T, 'SOC_Ah_Voltage_Current.csv');
disp('导出完成：SOC_Ah_Voltage_Current.csv');
