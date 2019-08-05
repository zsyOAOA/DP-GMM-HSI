clear;
clc;
close all;

% load data
load('./data/Terrain_std.mat');  % get Ori_H: 500

[M, N, B] = size(Ori_H);

r = 7;   % rank
initT_hdp = 5; initK_hdp = 30; iter_max_hap = 30;

% HDP
fprintf('\n')
disp('>>>>>>>>>>>>>>>>>>>>>>>>> Begin HDP >>>>>>>>>>>>>>>>>>>>>>>>>')
opts = set_opt_HSI('tol', 1e-4, 'initK', initK_hdp, 'initT', initT_hdp, 'itermax', iter_max_hap, 'init', 2, 'initR', r);
tic
[~, low_rank, varInfo] = hdp_denoise(Ori_H, opts);
toc
XHdp = low_rank.U * low_rank.V';
hdp_res = reshape(XHdp, [M, N, B]);


