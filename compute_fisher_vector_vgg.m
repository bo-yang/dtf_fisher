function [ fvt ] = compute_fisher_vector_vgg( desc, K , num_subsamp ,num_feats_per_video )
%COMPUTE_FISHER_VECTOR Compute Fisher vectors based on descriptors
%   Inputs:
%       desc - Nx1 dimention cell, collection of descriptors, each cell
%               element contains a channel of descriptor vectors(such as
%               HOG/HOF/MBHx/MBHy).
%       K - number of clusters(GMMs) to create
%       num_subsamp - number of subsampling vectors
%       num_feats_per_video - number of features per video
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

fvt=cell(length(desc)+1,1);
for i=1:length(desc)
    % Do PCA on train/test data to half-size original descriptors
    fprintf('Round %d: Performing PCA...\n', i);
    feats=desc{i};
    params.desc_dim = floor(size(feats,1)/2); % descriptor dimensionality after PCA
    low_proj = princomp(feats');
    low_proj = low_proj(:, 1:params.desc_dim)';
    if ~isempty(low_proj)
        % dimensionality reduction
        feats = low_proj * feats;
    end
    %feats=feats-repmat(mean(feats,1),size(feats,1),1);
    %sigma=feats*feats'/size(feats,2);
    %[U,~,~]=svd(sigma);
    %feats=U(:,1:params.desc_dim)'*feats;
    
    % Subsampling features for generaing GMM codebook
    idx = randperm(size(feats,2),num_subsamp);
    %idx = floor(linspace(1,size(faets,2),num_subsamp)); 
    v = feats(:,idx);   % random set of vectors for GMM estimation
    
    % Create codebook by GMM
    fprintf('Round %d: Creating codebook...\n',i);
    tic
    codebook=gmm_gen_codebook(v,params);
    toc
    
    fprintf('Round %d: Computing fish vector...\n',i);
    % Compute Fisher vectors
    for j=1:floor(size(feats,2)/num_feats_per_video)
        vt=feats(:,(j-1)*num_feats_per_video+1:j*num_feats_per_video);
        fvt{i} = [fvt{i} fisher_encode(vt,codebook,fisher_params)];
    end
    
    fvt{length(desc)+1}=[fvt{length(desc)+1};fvt{i}]; % the last one is combined descriptor
end

end

