function noise = multivariate_t_noise(R, nu, n_samples)
% Generate multivariate t-distributed noise
% Inputs:
%   R: target covariance matrix (d x d)
%   nu: degrees of freedom for t-distribution (must be > 2)
%   n_samples: number of noise samples
% Output:
%   noise: t-distributed noise matrix (n_samples x d)

% Check degrees of freedom
if nu <= 2
    error('Degrees of freedom nu must be greater than 2');
end

d = size(R, 1);  % Dimension

% Step 1: Adjust the covariance matrix for Gaussian distribution
adjusted_R = ((nu - 2) / nu) * R;

% Step 2: Generate multivariate Gaussian noise (zero mean)
gaussian_noise = mvnrnd(zeros(1, d), adjusted_R, n_samples);

% Step 3: Generate chi-squared random variables (with nu degrees of freedom)
chi2_samples = chi2rnd(nu, n_samples, 1);

% Step 4: Construct t-distributed noise
scaling_factors = sqrt(nu ./ chi2_samples);  % Scaling factors
noise = scaling_factors .* gaussian_noise;   % Scale sample-wise
end
