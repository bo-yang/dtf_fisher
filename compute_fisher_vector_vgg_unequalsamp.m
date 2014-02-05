function [ fvt ] = compute_fisher_vector_vgg( feats, K , num_subsamp ,feat_nums )
%COMPUTE_FISHER_VECTOR Compute Fisher vectors based on descriptors
%   Inputs:
%       desc - Nx1 dimention cell, collection of descriptors, each cell
%               element contains a channel of descriptor vectors(such as
%               HOG/HOF/MBHx/MBHy).
%       K - number of clusters(GMMs) to create
%       num_subsamp - number of subsampling vectors
%       feat_nums - vector, number of features per video clip
%
%   Outputs:
%       fvt - Nx1 dimention cell, collection of Fisher vectors

% To construct Fisher vector, fist to estimate the parameters of GMM.
% To estimate GMM params, subsampling vectors from DTF descriptors.

redo = 8;                  % number of times the initializing k-means is run (best clustering returned)
niter = 300;                % maximum number of iterations in k-means
verbose = 2;               % verbosity level
seed = 8;                  % 0: no seeding, values > 0 are used as seed
nt = 1;                    % to force multi-threading if not done by Matlab/octave
                           % check if multithreaded actived with nt=1 before changing this variable

eps=1e-5;

params.cluster_count=K;
params.maxcomps=params.cluster_count/4;
params.GMM_init= 'kmeans';


params.pnorm = single(2);    % L2 normalization, 0 to disable
fisher_params.alpha = single(0.5);  % power notmalization, 1 to disable
fisher_params.grad_weights = false; % soft BOW
fisher_params.grad_means = true;    % 1st order
fisher_params.grad_variances = true;    % 2nd order

params.subbin_norm_type = 'l2';
params.norm_type = 'l2';
params.post_norm_type = 'none';
params.pool_type = 'sum';
params.quad_divs = 2;
params.horiz_divs = 3;
params.kermap = 'hellinger';

fvt=zeros(size(feats,1)*K,numel(feat_nums));

% Preprocessing: L1 normalization & sqrt
feats=feats/norm(feats,1);
feats=sqrt(feats);

% Do PCA on train/test data to half-size original descriptors
fprintf('Round %d: Performing PCA...\n', i);
params.desc_dim = floor(size(feats,1)/2); % descriptor dimensionality after PCA
low_proj = princomp(feats');
low_proj = low_proj(:, 1:params.desc_dim)';
if ~isempty(low_proj)
    % dimensionality reduction
    feats = low_proj * feats;
end

% Subsampling features for generaing GMM codebook
idx = randperm(size(feats,2),num_subsamp);
%idx = floor(linspace(1,size(feats,2),num_subsamp));
v = feats(:,idx);   % random set of vectors for GMM estimation

% Create codebook by GMM
fprintf('Round %d: Creating codebook...\n',i);
tic
codebook=gmm_gen_codebook(v,params);
toc

fprintf('Round %d: Computing fish vector...\n',i);
% Compute Fisher vectors
start_idx=1;
for j=1:numel(feat_nums)
    vt=feats(:,start_idx:start_idx+feat_nums(j)-1);
    fvt(:,j) = fisher_encode(vt,codebook,fisher_params);
    start_idx=start_idx+feat_nums(j);
end

end

