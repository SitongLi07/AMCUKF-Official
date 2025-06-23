Adaptive Kernel Scale Unscented Kalman Filter (AMCUKF) Implementation
Nonlinear Numerical System （Section V.A）
File Structure
Core Algorithms
AUKF.m - Implementation of the proposed AMCUKF algorithm

EKF.m - Extended Kalman Filter

UKF_A1.m - Standard Unscented Kalman Filter

EUKF_iter.m - EMCUKF

FUKF_iter.m - FMCUKF

Main Execution
main.m - Primary script to reproduce all experimental results from Section V-A

main_diff_sigma.m - Generates kernel scale comparison plots (Fig. 3)


Noise Generation
cauchy_noise.m - Generates Cauchy-distributed noise

multivariate_t_noise.m - Generates multivariate t-distributed noise

Data Files
dis1.mat, dis2.mat, dis3.mat - Pre-generated disturbance profiles for experiments

Utility Functions
BaseZoom.m - Helper function for creating zoomed subplots


Getting Started
Basic Execution:

matlab
% Run all experiments from Section V-A
main();

% Generate kernel scale comparison plots (Fig. 3)
main_diff_sigma();
