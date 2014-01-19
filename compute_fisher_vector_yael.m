function [ fvt ] = compute_fisher_vector_yael( desc, K , num_subsamp ,num_feats_per_video )
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

fvt=cell(length(desc)+1,1);
for i=1:length(desc)
    % Do PCA on train/test data to half-size original descriptors
    desc{i}=removeDC(desc{i});    
    [V,~,~]=pca_mod(desc{i},eps);
    xw=V(1:floor(size(desc{i},1)/2),:)*desc{i};
    
    idx = randperm(size(xw,2),num_subsamp);
    %idx = floor(linspace(1,size(xw,2),num_subsamp)); 
    v = xw(:,idx);   % random set of vectors for GMM estimation
    
    tic
    [w, mu, sigma] = yael_gmm (single(v), K, 'redo', redo, 'niter', niter, 'seed', seed);
    toc
    
    % Compute a Fisher vector
    for j=1:floor(size(xw,2)/num_feats_per_video)
        vt=xw(:,(j-1)*num_feats_per_video+1:j*num_feats_per_video);
        fvt{i} = [fvt{i} yael_fisher(single(vt),w,mu,sigma,'sigma','nonorm')];
    end
    
    % power "normalization"
    fvt{i} = sign(fvt{i}) .* sqrt(abs(fvt{i}));
    
    % L2 normalization (may introduce NaN vectors)
    fvt{i} = yael_fvecs_normalize (single(fvt{i}));
    
    % replace NaN vectors with a large value that is far from everything else
    % For normalized vectors in high dimension, vector (0, ..., 0) is *close* to
    % many vectors.
    fvt{i}(find(isnan(fvt{i}))) = 123456;
    
    fvt{length(desc)+1}=[fvt{length(desc)+1};fvt{i}]; % the last one is combined descriptor
end

end

