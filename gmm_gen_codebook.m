function codebook = gmm_gen_codebook(feats,params)
%TRAIN Summary of this function goes here
%   Detailed explanation goes here

if (~exist('params','var'))
    params.cluster_count=256;
    params.numwords=4000000;
    params.samplesperclass=1500;
    params.numclases=60;
    params.maxcomps=params.cluster_count/4;
    params.desc_dim = 48; % descriptor dimensionality after PCA
    params.GMM_init= 'kmeans';
end

% -------------------------------------------------------------------------
%   Cluster codebook centres
% -------------------------------------------------------------------------

if isequal(params.GMM_init, 'kmeans')
    
    fprintf('Computing initial means using K-means...\n');

    % if maxcomps is below 1, then use exact kmeans, else use approximate
    % kmeans with maxcomps number of comparisons for distances
    if params.maxcomps < 1
        init_mean = vl_kmeans(feats, params.cluster_count, ...
            'verbose', 'algorithm', 'elkan');
    else
        init_mean = featpipem.lib.annkmeans(feats, params.cluster_count, ...
            'verbose', false, 'MaxNumComparisons', params.maxcomps, ...
            'MaxNumIterations', 150);
    end
    
    fprintf('Computing initial variances and coefficients...\n');

    % compute hard assignments
    kd_tree = vl_kdtreebuild(init_mean, 'numTrees', 3) ;
    assign = vl_kdtreequery(kd_tree, init_mean, feats);

    % mixing coefficients
    init_coef = single(vl_binsum(zeros(params.cluster_count, 1), 1, double(assign)));
    init_coef = init_coef / sum(init_coef);

    % variances
    init_var = zeros(size(feats, 1), params.cluster_count, 'single');

    for i = 1:params.cluster_count
        feats_cluster = feats(:, assign == i);
        init_var(:, i) = var(feats_cluster, 0, 2);
    end
    
elseif isequal(params.GMM_init, 'rand')
    init_mean = [];
    init_var = [];
    init_coef = [];
end

fprintf('Clustering features using GMM...\n');

% call FMM mex
gmm_params = struct;

if ~isempty(init_mean) && ~isempty(init_var) && ~isempty(init_coef)
    codebook = mexGmmTrainSP(single(feats), params.cluster_count, gmm_params, single(init_mean), single(init_var), single(init_coef));
else
    codebook = mexGmmTrainSP(single(feats), params.cluster_count, gmm_params);
end

fprintf('Done training codebook!\n');

end
