function x = cauchy_noise(x0, gamma, m, n)
% 生成 m×n 个柯西分布随机数
% x0: 位置参数 控制分布中心
% gamma: 尺度参数 越大，分布越“扁平”，尾部越厚
% m, n: 输出矩阵的行数与列数

    u = rand(m, n);  % 均匀分布 [0,1]
    x = x0 + gamma * tan(pi * (u - 0.5));
end
