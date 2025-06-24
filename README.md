markdown
# AMCUKF-Official

This repository contains MATLAB implementations for the paper:

**"Maximum Correntropy Unscented Kalman Filter with Unsupervised Adaptive Kernel Scale Selection"**

It includes three experimental scenarios:

---

## 1. Nonlinear Numerical System (Section V.A)

**Algorithms**  
- `AUKF.m` (Proposed AMCUKF)  
- `EKF.m`, `UKF_A1.m` (Baselines)  
- `EUKF_iter.m` (EMCUKF), `FUKF_iter.m` (FMCUKF)  

**Execution**
```matlab
main();             % Reproduce all Section V-A results  
main_diff_sigma();  % Generate Fig.3 kernel scale plots
```
Noise generation

cauchy_noise.m
multivariate_t_noise.m


## 2. Single-Target Tracking (Section V.B)

**Algorithms**  
Same as Scenario 1

**Execution**
```matlab
main();             % Reproduce all Section V-B results  
main_diff_sigma();  % Generate Fig.10 RMSE plots
```

## 3. Real-World SOC Estimation (Section V.C)

**Execution**
```matlab
UKF_results.m       % Reproduce all Section V-C results
```
For more details, please refer to the README.md files inside each scenario folder if available.
