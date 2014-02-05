function [ fvt ] = compute_fisher_vector_yael( feats, K, num_subsamp, feat_nums )
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

fvt=zeros(size(feats,1)*K,numel(feat_nums));


% Do PCA on train/test data to half-size original descriptors
feats=removeDC(feats);
[V,~,~]=pca_mod(feats,eps);
xw=V(1:floor(size(feats,1)/2),:)*feats;
%desc{i}=desc{i}-repmat(mean(desc{i},1),size(desc{i},1),1);
%sigma=desc{i}*desc{i}'/size(desc{i},2);
%[U,~,~]=svd(sigma);
%xw=U(:,1:floor(size(desc{i},1)/2))'*desc{i};

idx = randperm(size(xw,2),num_subsamp);
%idx = floor(linspace(1,size(xw,2),num_subsamp));
v = xw(:,idx);   % random set of vectors for GMM estimation

[w, mu, sigma] = yael_gmm (single(v), K, 'redo', redo, 'niter', niter, 'seed', seed);

% Compute a Fisher vector
start_idx=1;
for j=1:numel(feat_nums)
    vt=xw(:,start_idx:start_idx+feat_nums(j)-1);
    fvt(:,j) = yael_fisher(single(vt),w,mu,sigma,'sigma','nonorm');
    start_idx=start_idx+feat_nums(j);
end

% power "normalization"
fvt = sign(fvt) .* sqrt(abs(fvt));

% L2 normalization (may introduce NaN vectors)
fvt = yael_fvecs_normalize(single(fvt));

% replace NaN vectors with a large value that is far from everything else
% For normalized vectors in high dimension, vector (0, ..., 0) is *close* to
% many vectors.
fvt(find(isnan(fvt))) = 123456;

end

