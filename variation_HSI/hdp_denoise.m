function [param, uvinit, varInfo] = hdp_denoise(noiseData, opts, cleanData)
%Input:
%   noiseData: M X N X B tensor data.
%   opts: struct, returned by function set_opt.
%   cleanData: M x N x B original data
%Output:
%   param: MoG parameter, struct.
%   uvinit: Low Rank parameter. struct.

[M, N, B]      = size(noiseData);
d              = M*N;
Y              = reshape(noiseData, [d, B]);
param.sizedata = [M,N,B];
clear noiseData;

%Initialization
[param, uvinit, exValue] = hdp_initialize_HSI(Y, param, opts);

iter = 1;
currentTol = Inf;
bound = -Inf(opts.itermax+1, 1);
while (iter <= opts.itermax) && (currentTol > opts.tol)
    [param, exValue]= hdp_maximization_HSI(param, opts, exValue);
    [param, uvinit, exValue] = hdp_updateuv_HSI(Y, param, uvinit, exValue, opts);
    if opts.display
        LL       = uvinit.U * uvinit.V';
        err      = Y - LL;
        denoData = reshape(LL, [M N B]);
        if nargin == 3
            [mpsnr, mssim, ergas, ~, ~] = img_quality_HSI(denoData, cleanData);
            YClean = reshape(cleanData, [d, B]);
            GErr   = YClean - LL;
            fprintf('Iter:%2d, MPSNR=%.4f, MSSIM=%.4f, ERGAS=%08.4f, GEnorm=%08.4f,',...
                     iter, mpsnr, mssim, ergas, norm(GErr, 'fro'));
        else
            fprintf('Iter:%2d, Enorm=%.4f,',...
                     iter, norm(err,'fro'));
        end
    end
    [param, exValue] = hdp_expectation_HSI(Y, param, exValue);
    bound(iter) = low_bound_HSI(param, uvinit, opts, exValue);
    if iter > 1
        currentTol = abs(bound(iter)-bound(iter-1)) / abs(bound(iter));
    end
    if opts.display
        fprintf(' bound=%.4f, tol=%.4f...\n', bound(iter), currentTol);
    end
    iter = iter + 1;
end
varInfo = gauss_check(param);
param.bound = bound(1:iter);
end

function [param, uvinit, exValue] = hdp_initialize_HSI(Y, param, opts, cleanData)
% Input:
%   Y: d x B tensor data.
%   param: mog parameter
%   opts: hyperparameter
%   cleanData: M x N x B original data
% Output:
%   param: MoG parameter
%       param.rau: d x T x B tensor, responsibility matrix of local MoG.
%       param.pphi: T x K x B tensor, responsibility matrix of Global MoG.
%   uvinit: low rank matrix factorazation result
%       uvinit.U, uvinit.sigmaU;
%       uvinit.V, uvinit.sigmaV;

K      = opts.initK;
T      = opts.initT;
R      = opts.initR;
M      = param.sizedata(1);
N      = param.sizedata(2);
[d, B] = size(Y);
scale  = sqrt(sum(Y(:).^2)/(d*B));

switch opts.init
    case 1
        [U0, S, V0] = svd(Y, 'econ');
        U           = U0(:, 1:R) * S(1:R, 1:R)^0.6;
        Vt          = S(1:R, 1:R)^0.4 * V0(:, 1:R)';
        V           = Vt';
        if opts.display
            fprintf('SVD Initialization:');
        end
    case 2
        warning off;
        [A_hat, ~, ~] = inexact_alm_rpca(Y);
        [U, V]        = PCA(A_hat, R);
        warning on;
        if opts.display
            fprintf('RPCA Initialization:');
        end
    case 3
        [U, V]        = PCA(Y, R);
        if opts.display
            fprintf('PCA Initialization:');
        end
    case 4
        U   = rand(d, R) * sqrt(scale);
        V   = rand(B, R) * sqrt(scale);
        if opts.display
            fprintf('Random Initialization:');
        end
end
LL   = U * V';
err  = Y - LL;
if opts.display
    if nargin == 4
        value = zhibiao(cleanData, reshape(LL, [M, N, B]))
        fprintf('PSNR=%.4f,Enorm=%.4f...', value.mpsnr, norm(err, 'fro'));
    else
        fprintf('Enorm=%.4f...\n', norm(err, 'fro'));
    end
end

% Initialize U and V
uvinit.U          = U;
uvinit.V          = V;
uvinit.sigmaU     = repmat(eye(R)*scale, [1,1,d]);
uvinit.sigmaV     = repmat(eye(R)*scale, [1,1,B]);
exValue.errSquare = err.^2;

% Initialize q(C) and q(Z)
samplek = datasample(err(:), K, 'Replace', false); %K x 1 vector
pphi = zeros(T, K, B);
rau = zeros(d, T, B);
for bb = 1:B
    indexT = randsample(K, T);
    pphi(:,:,bb) = full(sparse(1:T, indexT, 1, T, K, T));
    rau(:,:,bb) = find_close_one(err(:, bb)', T, samplek(indexT)');
end
param.rau  = rau;
param.pphi = pphi;

% Initialize q(\alpha)
param.e       = opts.e0;
param.f       = opts.f0 * ones(1,B);
exValue.alpha = param.e ./ param.f;

% Initialize q(\gamma)
param.c        = opts.c0;
param.d        = opts.d0;
exValue.ggamma = param.c ./ param.d;

% Initialize q(\lambda)
r = 0.5*d + 0.5*B + opts.p0;
q = 0.5*diag(uvinit.U'*uvinit.U)+0.5*diag(sum(uvinit.sigmaU,3))...
            +0.5*diag(uvinit.V'*uvinit.V)+0.5*diag(sum(uvinit.sigmaV,3))...
            +opts.q0; % R x 1 vector
exValue.lambda = p ./ q;
end

function R = find_close_one(X, k, m)
% Input:
%   X: d x n matrix
%   m: d x k matrix
% Output:
%   R: n x K matrix, n is the number of data samples.
n = size(X, 2);
% m'*X: k x n matrix; dot(m,m,1): 1 x k vector
% label: 1 x n vector
[~, label] = max(bsxfun(@minus, m'*X, dot(m,m,1)'/2), [], 1);
[u, ~, label] = unique(label);
while k ~= length(u)
    idx = randsample(n, k);
    m = X(:, idx);
    [~, label] = max(bsxfun(@minus, m'*X, dot(m,m,1)'/2),[],1);
    [u, ~, label] = unique(label);
end
R = full(sparse(1:n,label,1,n,k,n));
end

function [param, exValue] = hdp_maximization_HSI(param, opts, exValue)

K         = size(param.pphi, 2);
[d, T, B] = size(param.rau);

%Update q(\pi')
numBT             = sum(permute(param.rau, [2, 3, 1]), 3); % T x B matrix
param.rrPi1       = numBT + 1;  % T x B matrix
param.rrPi2       = bsxfun(@minus, sum(numBT, 1), cumsum(numBT,1)) + repmat(exValue.alpha, [T, 1]);
exValue.piPie     = param.rrPi1 ./ (param.rrPi1+param.rrPi2); % T x B matrix
exValue.logPiPie  = psi(param.rrPi1) - psi(param.rrPi1+param.rrPi2); % T x B matrix
exValue.logPiPie1 = psi(param.rrPi2) - psi(param.rrPi1+param.rrPi2);% T x B matrix
temp              = [zeros(1, B); exValue.logPiPie1(1:end-1, :)]; % T x B matrix
exValue.logPi     = cumsum(temp, 1) + exValue.logPiPie; % T x B matrix
clear temp;

%Update q(beta')
ssBeta1             = reshape(permute(param.pphi, [2,1,3]), K, []) * ones(B*T, 1); % K x 1
param.ssBeta1       = ssBeta1 + 1; % K x 1 vector
param.ssBeta2       = sum(ssBeta1) - cumsum(ssBeta1) + exValue.ggamma; % K x 1
exValue.betaPie     = param.ssBeta1 ./ (param.ssBeta1+param.ssBeta2); % K x 1
exValue.logBetaPie  = psi(param.ssBeta1) - psi(param.ssBeta1+param.ssBeta2); % K x 1
exValue.logBetaPie1 = psi(param.ssBeta2) - psi(param.ssBeta1+param.ssBeta2); % K x 1
temp                = [0; exValue.logBetaPie1(1:end-1)];
exValue.logBeta     = cumsum(temp) + exValue.logBetaPie; % K x 1
clear ssBeta1;

%Update q(\gamma)
param.c          = K + opts.c0;
param.d          = opts.d0 - sum(exValue.logBetaPie1);
exValue.ggamma   = param.c / param.d;
exValue.logGamma = psi(param.c) - reallog(param.d);

%Update q(\alpha)
param.e          = T + opts.e0;
param.f          = opts.f0 - sum(exValue.logPiPie1, 1); % 1 x B vector
exValue.alpha    = param.e ./ param.f;
exValue.logAlpha = psi(param.e) - reallog(param.f);

%Update q(\xi)
numDataK = reshape(permute(param.pphi, [2,1,3]), K, [])...
            * reshape(permute(param.rau, [2,3,1]), [], d) * ones(d, 1); %K x 1 vector
param.a = 0.5*numDataK + opts.a0; %K x 1 vector
rauErrSquare = bsxfun(@times, permute(param.rau, [2,3,1]),...
                       permute(exValue.errSquare, [3,2,1])); % T x B x d matrix
rauPhiErrSquare =reshape(permute(param.pphi, [2,1,3]), K, []) * ...
                    reshape(rauErrSquare, [], d) * ones(d, 1); %K x 1 vector
param.b = opts.b0 + 0.5 * rauPhiErrSquare; %K x 1 vector
exValue.logXi = psi(param.a) - reallog(param.b); %K x 1 vector,E[ln(\xi_k)]
exValue.xi    = param.a ./ param.b;  %K x 1 vector
end

function [param, exValue] = hdp_expectation_HSI(Y, param, exValue)
% Update q(\rau) and q(C)

[d, B]  = size(Y);
[T, K, ~] = size(param.pphi);

%calculate E[ln N(y_ij-u_i*v_j|0, \xi_k)]
xiErr = bsxfun(@times, exValue.errSquare, reshape(exValue.xi, [1,1,K])); %d x B x K
exValue.logGauss = -0.5*(bsxfun(@minus, xiErr, reshape(exValue.logXi, [1,1,K])))...
             - 0.5 * log(2*pi); % d X B x K

%Update q(C)
pphi = zeros(T,K,B);
for jj = 1:B
    pphi(:,:,jj) = param.rau(:,:,jj)' * permute(exValue.logGauss(:,jj,:),[1,3,2]);
end
pphi = bsxfun(@plus, pphi, reshape(exValue.logBeta, [1,K,1]));
param.pphi = exp(bsxfun(@minus, pphi, logsumexp(pphi, 2)));
clear pphi;

%Update q(Z)
rau = zeros(d,T,B);
for jj = 1:B
    rau(:,:,jj) = permute(exValue.logGauss(:,jj,:), [1,3,2]) * param.pphi(:,:,jj)';
end
rau = bsxfun(@plus, rau, reshape(exValue.logPi, [1, T, B])); %d x T x B
rau = exp(bsxfun(@minus, rau, logsumexp(rau, 2)));
param.rau = rau;
end

function [param, uvinit, exValue] = hdp_updateuv_HSI(Y, param, uvinit, exValue, opts)
% Update U and V.

[d, B]    = size(Y);
R         = size(uvinit.U, 2);
[T, K, ~] = size(param.pphi);

phiXi = reshape(permute(param.pphi, [1,3,2]), [], K) * exValue.xi;
phiXi = reshape(phiXi, [T, B]); % T x B matrix
%Update U
VjVj = expectation_UV(uvinit.V, uvinit.sigmaV); % R x R x B
matVjVj = reshape(VjVj, R^2, B);
for ii=1:d
    rauii = permute(param.rau(ii,:,:), [2,3,1]); %T x B matrix
    rauPhiXk = sum(rauii .* phiXi, 1); %1 x B vector
    precisionii = reshape(matVjVj * rauPhiXk', R, R) + diag(exValue.lambda);
    uvinit.sigmaU(:,:,ii) = (precisionii)^-1;

    uvinit.U(ii,:) = (Y(ii,:).*rauPhiXk) * uvinit.V * uvinit.sigmaU(:,:,ii);
end
exValue.UrUr = diag(uvinit.U'*uvinit.U) + diag(sum(uvinit.sigmaU, 3));
clear VjVj rauii;

UiUi = expectation_UV(uvinit.U, uvinit.sigmaU); % R x R x d
matUiUi = reshape(UiUi, [], d);
%Updata V
for jj=1:B
    phiXk1 = param.pphi(:, :, jj) * exValue.xi; % T x 1 vector
    rauPhiXk = param.rau(:, :, jj) * phiXk1; % d x 1 vector
    precisionjj = reshape(matUiUi*rauPhiXk, R, R) + diag(exValue.lambda);
    uvinit.sigmaV(:,:,jj) = (precisionjj)^-1;

    uvinit.V(jj, :) = (rauPhiXk' .* Y(:,jj)') * uvinit.U * uvinit.sigmaV(:, :, jj);
end
exValue.VrVr = diag(uvinit.V'*uvinit.V) + diag(sum(uvinit.sigmaV, 3));
clear phiXk1 rauPhiXk;

VjVj = expectation_UV(uvinit.V, uvinit.sigmaV);
matVjVj = reshape(VjVj, R^2, B);
exValue.errSquare = Y.^2 + matUiUi' * matVjVj - 2*Y.*(uvinit.U * uvinit.V');
clear matUiUi matVjVj;

%Update q(\lambda)
param.p = 0.5*d + 0.5*B + opts.p0;
param.q = 0.5*exValue.UrUr + 0.5 * exValue.VrVr + opts.q0;
exValue.lambda = param.p ./ param.q;
exValue.logLambda = psi(param.p) - reallog(param.q);

if opts.cropR && R > opts.MinR
    index = find(exValue.lambda<=opts.thresUV);
    if length(index) < opts.MinR
        [~, index0] = sort(exValue.lambda, 'ascend');
        validRIndex = index0(1:opts.MinR);
    else
        validRIndex = index;
    end
    uvinit.U = uvinit.U(:, validRIndex);
    uvinit.V = uvinit.V(:, validRIndex);
    uvinit.sigmaU = uvinit.sigmaU(validRIndex, validRIndex, :);
    uvinit.sigmaV = uvinit.sigmaV(validRIndex, validRIndex, :);
    param.q = param.q(validRIndex);
    param.exLambda = param.exLambda(validRIndex);
end
end

function bound = low_bound_HSI(param, uvinit, opts, exValue)
% Calculate the variation lower bound.

[T, K, B] = size(param.pphi);
R         = size(uvinit.U, 2);
M         = param.sizedata(1);
N         = param.sizedata(2);
d         = M * N;

bound = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%calculate E[q(W)log(p(X,W))]%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%E[ln(p(Y|U,V,\theta,C, Z))]
phiErr = zeros(d,T,B);
for jj = 1:B
    phiErr(:,:,jj) = permute(exValue.logGauss(:,jj,:), [1,3,2])...
                             * param.pphi(:, :, jj)';
end
rauPhiGauss = phiErr .* param.rau; %d x T x B
bound = bound + sum(rauPhiGauss(:));
clear rauPhiGauss phiErr;

%E[ln(p(U|\lambda))] and E[ln(p(V|\lambda))]
bound = bound + 0.5*sum(d*exValue.logLambda - d*log(2*pi)...
              - exValue.lambda.*exValue.UrUr);
bound = bound + 0.5*sum(B*exValue.logLambda - B*log(2*pi)...
              - exValue.lambda.*exValue.VrVr);

%E[ln(p(\lambda))]
bound = bound + sum(opts.p0*log(opts.q0) - gammaln(opts.p0)...
              + (opts.p0-1)*exValue.logLambda - opts.q0*exValue.lambda);

%E[p(\theta)]
bound = bound + sum(opts.a0*log(opts.b0) - gammaln(opts.a0) ...
              + (opts.a0-1)*exValue.logXi - opts.b0*exValue.xi);

%E[ln(p(C))]
phiExLogBeta = reshape(permute(param.pphi, [1,3,2]), [], K) * exValue.logBeta;
bound = bound + sum(phiExLogBeta);
clear phiExLogBeta;

%E[ln(p(Beta^{'}))]
bound = bound + sum(exValue.logGamma + (exValue.ggamma-1)*exValue.logBetaPie1);

%E[ln(p(gamma))]
bound = bound + opts.c0*log(opts.d0) - gammaln(opts.c0)...
              + (opts.c0-1)*exValue.logGamma - opts.d0*exValue.ggamma;

%E[ln(p(Z))]
rauExLogPi = bsxfun(@times, param.rau, reshape(exValue.logPi, [1,T,B]));
bound = bound + sum(rauExLogPi(:));
clear rauExLogPi;

%E[ln(p(pi^{'}))]
alphaPi = bsxfun(@times, exValue.logPiPie1, exValue.alpha-1);  %T x B matrix
bound = bound + sum(sum(bsxfun(@plus, alphaPi, exValue.logAlpha)));
clear exLogPiPie1;

%E[ln(p(alpha))]
bound = bound + sum(opts.e0*log(opts.f0) - gammaln(opts.e0)...
              + (opts.e0-1)*exValue.logAlpha - opts.f0*exValue.alpha);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lcalculate E[q(W)log(q(W))]%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%E[lnq(Z)]
logRau = log(param.rau);
logRau(isinf(logRau)) = 0;
rauLogRau = param.rau .* logRau;
bound = bound - sum(rauLogRau(:));
clear rauLogRau logRau;

%E[lnq(C)]
logPhi = log(param.pphi);
logPhi(isinf(logPhi)) = 0;
phiLogPhi = param.pphi .* logPhi;
bound = bound - sum(phiLogPhi(:));
clear phiLogPhi logPhi;

%E[lnq(pi^{'})]
bLog = betaln(param.rrPi1, param.rrPi2); % T x B matrix
entropyPi = bLog - (param.rrPi1-1) .* psi(param.rrPi1)...
                 - (param.rrPi2-1) .* psi(param.rrPi2)...
                 + (param.rrPi1+param.rrPi2-2) .* psi(param.rrPi1+param.rrPi2);
bound = bound + sum(entropyPi(:));

%E[lnq((Beta^{'}))]
bLogss = betaln(param.ssBeta1, param.ssBeta2);
entropyBetaPie = bLogss - (param.ssBeta1-1) .* psi(param.ssBeta1)...
                        - (param.ssBeta2-1) .* psi(param.ssBeta2)...
                        + (param.ssBeta1+param.ssBeta2-2)...
                                .* psi(param.ssBeta1+param.ssBeta2);
bound = bound + sum(entropyBetaPie(:));

%E[lnq(gamma)]
bound = bound + param.c - reallog(param.d) + gammaln(param.c)...
              + (1-param.c) * psi(param.c);

%E[ln(q(\alpha))]
bound = bound + sum(param.e - reallog(param.f) + gammaln(param.e)...
              + (1-param.e).*psi(param.e));

%E[q(\theta)]
bound = bound + sum(param.a - reallog(param.b) + gammaln(param.a)...
              + (1-param.a) .* psi(param.a));

%E[ln(q(U))]
for ii = 1:d
    cholUi = chol(uvinit.sigmaU(:, :, ii));
    logDetUi = sum(reallog(diag(cholUi)));
    bound = bound + 0.5*R*(1+log(2*pi)) + logDetUi;
end

%E[ln(q(V))]
for jj=1:B
    cholVj = chol(uvinit.sigmaV(:, :, jj));
    logDetVj = sum(reallog(diag(cholVj)));
    bound = bound + 0.5*R*(1+log(2*pi)) + logDetVj;
end

%E[ln(q(\lambda))]
bound = bound + sum(param.p - reallog(param.q) + gammaln(param.p) ...
              + (1-param.p) .* psi(param.p));
end

function varInfo = gauss_check(mogParam)
% Determine the number of gauss component.
% Input:
%   mogParam: struct, returned by function hdp_multi_view
%       mogParam.rau: d x T x B tensor, d =M*N
%       mogParam.phi: T x K x B tensor
% Output:
%   globalVar: the gauss variance of global MoG
%   groupVar: the gauss variance of each group MoG

[d, T, B] = size(mogParam.rau);
variance    = mogParam.b ./ (mogParam.a-1);
K            = length(variance);
rau          = mogParam.rau; % d x T x B tensor
phi          = mogParam.pphi; % T x K x B tensor

try
    ggamma = mtimesx(rau, phi); % d x K x B tensor
catch
    ggamma = zeros(d, K, B);
    for jj = 1:B
        ggamma(:, :, jj) = rau(:, :, jj) * phi(:, :, jj);
    end
end
ggamma = bsxfun(@rdivide, ggamma, sum(ggamma, 2)); % d x K x B tensor

groupVar  = cell(B, 1);
indexVar  = cell(B, 1);
mask      = zeros(d, B);
for jj = 1:B
    [~, indjj0]    = max(ggamma(:, :, jj), [], 2);
    mask(:, jj) = indjj0;
    indjj          = unique(indjj0);
    [~, ind_sort_var] = sort(variance(indjj));
    indjj = indjj(ind_sort_var);
    groupVar{jj}   = variance(indjj);
    indexVar{jj}   = indjj;
end
varInfo.globalVar = variance;
varInfo.mask      = mask;
varInfo.groupVar  = groupVar;
varInfo.indexVar  = indexVar;
end

function exXiX = expectation_UV(X, sigmaX)
%Input:
%   X: d x R matrix
%Output:
%   exXiX: R x R x d, E[u_{i \cdot}^T * u_{i \cdot}]
[d, R] = size(X);
XasTensor1 = reshape(X', [R,1,d]); %R x 1 x d
XasTensor2 = reshape(X', [1,R,d]);
exXiX = bsxfun(@times, XasTensor1, XasTensor2) + sigmaX; % R x R x d
end

