clear;
clc;

data_path = './data/pure_DCmall.mat';
load(data_path);   % get variable Ori_H

[M, N, B] = size(Ori_H);

r = 7;   % rank
seed = 10000;

case_num = 2;
band_deadline = [];
band_stripe = [];
switch case_num
    case 1
        noi_H = noise_case1(Ori_H, 0.0025, seed);
    case 2
        noi_H = noise_case2(Ori_H, seed);
    case 3
        [noi_H, band_deadline]= noise_case3(Ori_H, 40, seed, 1);
    case 4
        [noi_H, band_stripe]= noise_case4(Ori_H, 30, seed, 1);   
    case 5
        [noi_H]= noise_case5(Ori_H, seed);   
    case 6
        [noi_H, band_deadline, band_stripe]= noise_case6(Ori_H, 20, 20, seed);
end

noi_H_mat = reshape(noi_H, [], B);

% save result
%        1.noise    2.HDP
% mpsnr  
% mssim  
% egras  
quality = zeros(3, 2);
[quality(1,1), quality(2,1), quality(3,1), noi_psnr, noi_ssim] = img_quality_HSI(noi_H, Ori_H);

% HDP
fprintf('Case %d\n', case_num)
disp('>>>>>>>>>>>>>>>>>>>>>>>>> Begin HDP >>>>>>>>>>>>>>>>>>>>>>>>>')
opts = set_opt_HSI('tol', 1e-4, 'initK', 40, 'initT', 5, 'itermax', 30, 'init', 2, 'initR', r);
tic
[~, low_rank, varInfo] = hdp_denoise(noi_H, opts, Ori_H);
toc
XHdp = low_rank.U * low_rank.V';
hdp_res = reshape(XHdp, [M, N, B]);
[quality(1, 2), quality(2, 2), quality(3, 2), hdp_psnr, hdp_ssim] = img_quality_HSI(hdp_res, Ori_H);

fprintf('\nNoisy:    MPSNR=%4.2f, MSSIM=%5.4f, EGRAS=%06.2f\n', quality(1,1), quality(2,1), quality(3,1))
fprintf('Denoised: MPSNR=%4.2f, MSSIM=%5.4f, EGRAS=%06.2f\n', quality(1,2), quality(2,2), quality(3,2))

