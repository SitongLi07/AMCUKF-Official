
% AMCUKF

function [SOC_est, V_est, P_trace, state_history,iter_num] = simulate_ESC_AUKF(X, current, voltage_meas, soc0)
   
    alpha = 1e-3; kappa = 0; beta = 2;
    L = 3;  
    lambda = alpha^2 * (L + kappa) - L;
    gamma = sqrt(L + lambda);

  
    Wm = [lambda / (L + lambda), repmat(1 / (2 * (L + lambda)), 1, 2*L)];
    Wc = Wm;
    Wc(1) = Wc(1) + (1 - alpha^2 + beta);

   
    deltaT = 60 / 3600;
    eta = 1;
    N = length(current);
    x = [soc0; 0; 0];
    P = eye(3) * 1e-4;

    Q = diag([1e-7, 1e-6, 1e-7]);
    R = 1e-4;

    SOC_est = zeros(N,1);
    V_est = zeros(N,1);
    P_trace = zeros(N,1);
    state_history = zeros(N, 3);
    s_prev = sign(current(1));

    aalpha0 = 3; bbeta0 = 3; rhoA = 0.9; rhoB = 0.9;
    aalpha = aalpha0; bbeta = bbeta0;

    for k = 1:N
        i = -current(k);

      
        S = chol(P, 'lower');
        X_sigma = [x, x + gamma*S, x - gamma*S];

   
        X_pred = zeros(L, 2*L+1);
        for j = 1:2*L+1
            X_pred(:, j) = esc_state_transition(X_sigma(:, j), i, deltaT, eta, X, s_prev);
        end

        x_pred = X_pred * Wm';
        P_minus = Q;
        for j = 1:2*L+1
            dx = X_pred(:, j) - x_pred;
            P_minus = P_minus + Wc(j) * (dx * dx');
        end

     
        Z_sigma = zeros(1, 2*L+1);
        for j = 1:2*L+1
            Z_sigma(j) = esc_measurement_model(X_pred(:, j), i, X, s_prev);
        end

        z_pred = Z_sigma * Wm';

        Pxz = zeros(L, 1);
        for j = 1:2*L+1
            Pxz = Pxz + Wc(j) * (X_pred(:, j) - x_pred) * (Z_sigma(j) - z_pred);
        end

        x = x_pred;
        iter = 0;
        aalpha_minus = rhoA * aalpha; 
        bbeta_minus = rhoB * bbeta;
        aalpha = 1/2 + aalpha_minus; 
        bbeta = bbeta_minus;
        while(iter < 5)
            xtemp = x;
            tau = bbeta / (aalpha - 1);
            Pz = tau * R;
            for j = 1:2*L+1
                dz = Z_sigma(j) - z_pred;
                Pz = Pz + Wc(j) * dz^2;
            end

            K = Pxz / Pz;
            y = voltage_meas(k) - z_pred;
            x = x_pred + K * y;
            P = P_minus - Pxz / Pz * Pxz';


            z_kk = esc_measurement_model(x, i, X, s_prev);

            y_err = voltage_meas(k) - z_kk;
            bbeta = bbeta_minus + 0.5*(y_err)' * inv(R) * (y_err);
            iter = iter + 1;
            if (abs(x - xtemp) / abs(xtemp) <= 1e-3)
                break;
            end
        end

      
        SOC_est(k) = x(1);
        V_est(k) = z_pred;
        P_trace(k) = trace(P);
        state_history(k, :) = x';
        iter_num(k) = iter;
        tau_t(k) = bbeta / (aalpha -1);

        if abs(i) > 0
            s_prev = sign(i);
        end
    end
end


function x_next = esc_state_transition(x, i, deltaT, eta, X, s_prev)
    M0 = X(1); M = X(2); gamma = X(3);
    R0 = X(4); R1 = X(5); Cn = X(6);

    ARC = exp(-deltaT / (R1 * Cn));
    BRC = 1 - ARC;
    AH = exp(-abs(-eta * gamma * deltaT * i / Cn));
    BH = 1 - AH;

    SOC_k1 = x(1) - (deltaT / Cn) * eta * i;
    iR_k1  = ARC * x(2) + BRC * i;
    h_k1   = AH * x(3) - BH * s_prev;

    x_next = [SOC_k1; iR_k1; h_k1];
end


function z = esc_measurement_model(x, i, X, s_prev)
    M0 = X(1); M = X(2);
    R0 = X(4); R1 = X(5);

    soc = x(1); iR = x(2); h = x(3);
    OCV = ocv_from_soc(soc);
    z = OCV + M0 * s_prev + M * h + R1 * iR + R0 * i;
end
