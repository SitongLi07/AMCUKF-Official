Adaptive Kernel Scale Unscented Kalman Filter (AMCUKF) Implementation
Real-World SOC Estimation （Section V.A）

Raw and Corrected Data:

10-28-18_14.16 551_Charge2.mat — Original charging data ---------------------- Data 1

SOC_Ah_Voltage_Current.csv — Zero-offset corrected SOC, Ah, voltage, and current data -------- Data 2

Fitting:

export_soc_ah_voltage_current.m — Loads Data 1, generates Data 2; Data 2 is the result of zero-offset correction applied to Data 1

plot_soc_comparison.m — Plots the SOC after zero correction

fit_ocv_soc_poly6.m — Performs 6th-order polynomial fitting, displays the fitted curve, and outputs OCV(SOC) polynomial coefficients

Parameter Identification:

simulate_ESC_model.m — Simulates state updates and voltage estimation using the ESC model

ocv_from_soc.m — Explicit 6th-order polynomial expression for OCV(SOC)

esc_objective_function.m — Minimizes the mean squared error (MSE) between model-estimated voltage and measured voltage

identify_ESC_parameters.m — Main script to obtain optimal ESC model parameters

Estimation:

simulate_ESC_UKF.m — UKF algorithm, encapsulated as a function

simulate_ESC_EUKF.m — EMCUKF algorithm, encapsulated as a function

simulate_ESC_FUKF.m — FMCUKF algorithm, encapsulated as a function

simulate_ESC_AUKF.m — AMCUKF algorithm, encapsulated as a function

UKF_results.m — Main script to produce experimental results in Section V.C

Utility Functions:

BaseZoom.m — Helper function for creating zoomed-in subplots
