function Y = noise_case1(X, variance, seed)
% IID Gaussian noise.
% Input:
%   X: M x N X B tensor
%   variance: variance of Gaussian noise, default 0.02
%   seed: seed of random number, default: 10000
% Output:
%   Y: M x N x B noise tensor

if nargin == 1
    variance = 0.0025;
    seed = 10000;
elseif nargin == 2
    seed = 10000;
end

[M, N, B] = size(X);
rng('default'); rng(seed);

Y = zeros(M, N, B);
for bb = 1:B
    Y(:, :, bb) = X(:, :, bb) + sqrt(variance) * randn(M, N);
end

end
