function Y = noise_case5(X, seed)
% Add impluse noise to HSIs
% Input:
%   X: M x N x B tensor.
%   seed: seed of random number, default 10004
% Output:
%   Y: M x N x B noise tensor

if nargin == 1
    seed = 10004;
end

B = size(X, 3);
Y = noise_case2(X, seed);

rng('default'); rng(seed);
intensity = randsample(0:0.01:0.15, B, 'true');

for jj = 1:B
    Y(:,:,jj) = imnoise(Y(:, :, jj), 'salt & pepper', intensity(jj));
end

end

