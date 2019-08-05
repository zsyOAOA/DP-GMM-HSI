function psnr = psnr_img(noise, clean)
% Calculate the PSNR value between the noisy image and the clean image.
% Input:
%   noise: h x w x 3 RGB image or h x w gray image (noise), uint8 range [1, 255]
%   clean: h x w x 3 RGB image or h x w gray image (clean), uint8 range [1, 255]
% Output:
%   psrn: PSNR value

if any(isnan(noise))
    error('The noise image contains NAN value');
end
if any(isnan(clean))
    error('The clean image contains NAN value');
end

num_dim = ndims(noise);
h       = size(noise, 1);
w       = size(noise, 2);

if num_dim == 3
    diff_square = mean((noise-clean).^2, 3);
else
    diff_square = (noise - clean).^2;
end

mse = sum(diff_square(:)) / (h*w);

psnr = 20 * log10(255) - 10 * log10(mse);

end

