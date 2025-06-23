Adaptive Kernel Scale Unscented Kalman Filter (AMCUKF) Implementation
Single-Target Tracking Example （Section V.B）
File Structure
Core Algorithms
AUKF.m - Implementation of the proposed AMCUKF algorithm

EKF.m - Extended Kalman Filter

UKF_A1.m - Standard Unscented Kalman Filter

EUKF_iter.m - EMCUKF

FUKF_iter.m - FMCUKF

Main Execution
main.m - Primary script to reproduce all experimental results from Section V-A

main_diff_sigma.m - Generates The statistical average of the position RMSE over time of the FMCUKF with different kernel scale values (Fig.10)


Noise Generation
cauchy_noise.m - Generates Cauchy-distributed noise

multivariate_t_noise.m - Generates multivariate t-distributed noise

Data Files
dis.mat - Pre-generated disturbance profiles for experiments

Utility Functions
BaseZoom.m - Helper function for creating zoomed subplots


Getting Started
Basic Execution:

matlab
% Run all experiments from Section V-A
main();

% Generates The statistical average of the position RMSE over time of the FMCUKF with different kernel scale values (Fig.10)
main_diff_sigma();