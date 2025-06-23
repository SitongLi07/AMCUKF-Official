
% ESC model simulator

function [V_est, state_history] = simulate_ESC_model(X, current, soc0)
    % Input parameters:
    % X = [M0, M, gamma, R0, R1, Cn]';
    % current: current sequence (Nx1), unit: A
    % soc0: initial SOC value (scalar)
    % Output:
    % V_est: simulated estimated voltage sequence (Nx1)
    % state_history: state at each time step [SOC, iR, h]

   
    M0 = X(1);
    M = X(2);
    gamma = X(3);
    R0 = X(4);
    R1 = X(5);
    Cn = X(6);

    N = length(current);         % Data length
    deltaT = 60 / 3600;          % Sampling period (in hours)
    eta = 1;                     % Coulombic efficiency, assumed constant

    % 初始状态 [SOC, iR, h]
    state = [soc0; 0; 0];
    state_history = zeros(N, 3);
    V_est = zeros(N, 1);
    s_prev = sign(current(1));

    for k = 1:N
        i = -current(k);  

      
        ARC = exp(-deltaT / (R1 * Cn));
        BRC = 1 - ARC;

        AH = exp(-abs(-eta * gamma * deltaT * i / Cn));
        BH = 1 - AH;

        
        SOC_k1 = state(1) - (deltaT / Cn) * eta * i;
        iR_k1  = ARC * state(2) + BRC * i;
        h_k1   = AH * state(3) - BH * s_prev;

    
        state = [SOC_k1; iR_k1; h_k1];
        state_history(k, :) = state';

     
        hyst = M0 * s_prev + M * h_k1;

        
        OCV = ocv_from_soc(SOC_k1);  
        V_est(k) = OCV + hyst + R1 * iR_k1 + R0 * i;

      
        if abs(i) > 0
            s_prev = sign(i);
        end
    end
end
