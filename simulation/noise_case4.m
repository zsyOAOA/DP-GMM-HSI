function [Y, band_stripe]= noise_case4(X, num_bands, seed, gauss)
% Add stripe noise to HSIs.
% Input:
%   X: M x N x B tensor
%   num_bands: number of bands that are added stripe noise
%   seed: seed of random number, default 10002
%   gauss: 0 or 1
% Output:
%   Y: M x N x B tensor

if nargin == 2
    seed = 10002;
    gauss = 1;
end
if nargin == 3
    gauss = 1;
end

[M, ~, B] = size(X);
if gauss == 1
    Y = noise_case2(X, seed);
else
    Y = X;
end
max_stripe_per_band = 40;

rng('default'); rng(seed);

band_stripe = randperm(B, num_bands)';
num_stripe = randi([15,max_stripe_per_band], num_bands, 1);
width_stripe = randi([1,3], num_bands, max_stripe_per_band);

for jj = 1:num_bands
    try
        location = randsample(1:4:M-4, num_stripe(jj));
    catch
        location = randsample(1:4:M-4, num_stripe(jj), 'true');
        location = unique(location);
    end
    for ii = 1:length(location)
        loc_ii = location(ii);
        t = rand(width_stripe(jj,ii), 1) * 0.5 - 0.15;
        Y(loc_ii:loc_ii+width_stripe(jj,ii)-1, :, band_stripe(jj)) = bsxfun(@minus, Y(loc_ii:loc_ii+width_stripe(jj,ii)-1, :, band_stripe(jj)), t);
    end
end

end

