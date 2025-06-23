
% File: esc_objective_function.m

function J = esc_objective_function(X, current, voltage_meas, soc0)
% Input:
%   X - Parameter vector [M0, M, gamma, R0, R1, Cn]
%   current - Current vector (in Amperes)
%   voltage_meas - Measured voltage vector
%   soc0 - Initial SOC

% Estimate voltage using ESC model simulator
    try
        V_est = simulate_ESC_model(X, current, soc0);
    catch
        J = 1e6; 
        return;
    end

    if length(V_est) ~= length(voltage_meas)
        J = 1e6;
    else
        error_vec = voltage_meas - V_est;
        J = mean(error_vec .^ 2);
    end
end
