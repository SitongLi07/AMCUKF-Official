function noise = multivariate_t_noise(R, nu, n_samples)
% 生成多变量t分布噪声
% 输入:
%   R: 目标协方差矩阵 (d x d)
%   nu: t分布自由度 (必须 > 2)
%   n_samples: 噪声样本数
% 输出:
%   noise: t分布噪声矩阵 (n_samples x d)

% 检查自由度
if nu <= 2
    error('自由度 nu 必须大于 2');
end

d = size(R, 1);  % 维度

% 步骤1: 调整高斯分布的协方差矩阵
adjusted_R = ((nu - 2) / nu) * R;

% 步骤2: 生成多元高斯噪声 (均值为0)
gaussian_noise = mvnrnd(zeros(1, d), adjusted_R, n_samples);

% 步骤3: 生成卡方随机变量 (nu自由度)
chi2_samples = chi2rnd(nu, n_samples, 1);

% 步骤4: 组合成t分布噪声
scaling_factors = sqrt(nu ./ chi2_samples);  % 缩放因子
noise = scaling_factors .* gaussian_noise;   % 按样本缩放
end