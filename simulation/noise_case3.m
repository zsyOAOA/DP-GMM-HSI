function [Y, band_deadline]= noise_case3(X, num_bands, seed, gauss)
% Add deadlines to HSIs.
% Input:
%   X: M x N x B tensor
%   numbands: the num of bands that are added noise
%   seed: seed of random number, default 10001
%   gauss: 0 or 1
% Output:
%   Y: M x N x B tensor

if nargin == 1
    num_bands = 40;
    seed = 10001;
    gauss = 1;
end
if nargin == 2
    seed = 10001;
    gauss = 1;
end
if nargin == 3
    gauss = 1;
end

[~, N, B] = size(X);
if gauss == 1
    Y = noise_case2(X, seed);
else
    Y = X;
end
max_deadline_per_band = 15;

rng('default'); rng(seed);

band_deadline = randperm(B, num_bands)';
num_deadline = randi([5,max_deadline_per_band], num_bands, 1);
width_deadline = randi([1,2], num_bands, max_deadline_per_band);

for jj = 1:num_bands
    try
        location = randsample(1:4:N-4, num_deadline(jj), 'false');
    catch
        location = randsample(1:4:N-4, num_deadline(jj), 'true');
        location = unique(location);
    end
    for ii = 1:length(location)
        loc_ii = location(ii);
        Y(:, loc_ii:loc_ii+width_deadline(jj, ii)-1, band_deadline(jj)) = 0;
    end
end

end

