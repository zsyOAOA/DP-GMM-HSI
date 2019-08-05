function opt = set_opt_HSI(varargin)
%set_opt.m set parameters of the HDP MoG Model,return a structure containing complete parameter setting
% usage:
%      opt = set_opt('Param1','StringValue1','Param2',NumValue2)
% Output:
%      opt: struct of hyperparameter.

opt = struct('a0', 1e-6,...
             'b0', 1e-6,...      %\xi \sim Gamma(a0, b0)
             'c0', 1e-6,...
             'd0', 0.25*1e-6,... %\gamma \sim Gamma(c0, d0), scaling parameter of Global Dirichlet Process
             'e0', 1e-6,...
             'f0', 0.5*1e-6,...  %\alpha \sim Gamma(e0, f0), scaling parameter of local Dirichlet Process
             'p0', 1e-3,...
             'q0', 1e-3,...      %\lambda \sim Gamma(p0, q0)
             'initK', 40,...     %Initialization number of Global MoG components
             'initT', 5,...      %Initialization number of local MoG components
             'cropR', 0,...
             'initR', 5,...
             'MinR',5,...
             'thresUV',1e4,...
             'tol', 1e-6,...
             'itermax', 20,...
             'display', 1,...
             'init',2);          % Initialization of U and V: 1:svd, 2:RPCA, 3:PCA, 4:random
numparameters = nargin;
if numparameters == 0
    return
end
if rem(numparameters, 2) ~= 0
    error('The paarmeters must be field-value pairs')
end
if rem(numparameters, 2) == 0 && numparameters > 1
    for i = 1:2:numparameters
        opt.(varargin{i}) = varargin{i+1};
    end
end

