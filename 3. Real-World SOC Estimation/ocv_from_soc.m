
% File: ocv_from_soc.m

function OCV = ocv_from_soc(soc)
    % Input:
    % soc: SOC values (can be a vector)
    % Output:
    % OCV: Voltage values calculated using the fitted polynomial

    % 6th-order polynomial coefficients (from highest to lowest degree)
    %coeffs = [12.2152, -26.9770, 12.1631, 9.3358, -8.9846, 3.1420, 3.3051];
    coeffs = [10.6258,-17.6302, -4.2855, 21.5811, -13.1342, 3.7775, 3.2818];
       
    % 计算多项式值
    OCV = polyval(coeffs, soc);
end
