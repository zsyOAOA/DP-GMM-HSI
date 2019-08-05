function [mpsnr, mssim, ergas, psnr_val, ssim_val] = img_quality_HSI(noise_H, ref_H)
% image quality assessment of SSIM, PSNR, ERGAS for HSIs data.
% Input:
%   noise_H: M x N x B HSI data (noise), range [0,1]
%   ref_H: M x N x B HSI data (reference), range [0,1]
% Output:
%   mssim: scalar, MSSIM value
%   mpsnr: scalar, MPSNR value
%   egras: scalar, SRGAS value
%   ssim_val: B x 1 vector, SSIM value
%   psnr_val: B x 1 vector, PSNR value

noise_H = noise_H * 255;
ref_H = ref_H * 255;

ergas = ErrRelGlobAdimSyn(noise_H, ref_H);

B = size(noise_H, 3);
ssim_val = zeros(B, 1);
psnr_val = zeros(B, 1);
for jj = 1:B
    [ssim_val(jj), ~] = ssim_index(noise_H(:, :, jj), ref_H(:, :, jj));
    psnr_val(jj) = psnr_img(noise_H(:, :, jj), ref_H(:, :, jj));
end
mssim = mean(ssim_val);
mpsnr = mean(psnr_val);

end

