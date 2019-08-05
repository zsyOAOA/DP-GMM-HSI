function [Y, band_deadline, band_stripe]= noise_case6(X, num_bands_stripe, num_bands_deadline, seed)
% Add stripe, deadline, impulse noise to HSIs.
% Input:
%   X: M x N x B tensor
%   num_bands_stripe: number of bands that are added stripe noise
%   num_bands_deadline: number of bands that are added deadlines
%   seed: seed of random number
% Output:
%   Y: M x N noise tensor

if nargin == 3
    seed = 10005;
end

Y_5 = noise_case5(X, seed);

[Y_5_3, band_deadline]= noise_case3(Y_5, num_bands_deadline, seed+1, 0);

[Y, band_stripe]= noise_case4(Y_5_3, num_bands_stripe, seed+2, 0);

end



