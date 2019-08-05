function [Y, variance]= noise_case2(X, seed)
% Non-i.i.d. zero-mean Gaussion noise, variance of each band is select from 0.01:0.003:0.04.
% Input:
%   X: M x N x B tensor
%   seed: seed of random number, default 10000
% Output:
%   Y: M x N x B tensor

[M, N, B] = size(X);

if nargin == 1
    seed = 10000;
end

rng('default'); rng(seed);

var_space = 0.001:0.001:0.01;
variance = randsample(var_space, B, true);

Y = zeros(M, N, B);
for jj = 1:B
    Y(:, :, jj) = X(:, :, jj) + sqrt(variance(jj)) * randn(M, N);
end

end
