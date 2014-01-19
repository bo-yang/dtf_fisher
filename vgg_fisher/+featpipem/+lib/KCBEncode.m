function encoding = KCBEncode(feats, ...
    dictwords, K, sigma, kdtree, maxComparisons, kcbType, outputFullMat)
%KCBENCODE Summary of this function goes here
%   Detailed explanation goes here

    % START validate input parameters -------------------------------------
%     default('K', 5);
%     default('sigma', 1e-4);
    % only used if a kdtree is specified
%     default('maxComparisons', size(dictwords,2));
    % possible values: 'unc', 'pla', 'kcb'
    % see van Gemert et al. ECCV 2008 for details
    % 'inveuc' simply weights by the inverse of euclidean distance, then l1
    % normalizes
%     default('kcbType', 'unc');
%     default('outputFullMat', true);
    % END validate input parameters ---------------------------------------
    
    % if using the 'plausibility' method, only 1NN is required
    if strcmp(kcbType, 'pla')
        K = 1;
    end
    
    num_words = size(dictwords, 2);
    num_feats = size(feats, 2);
    
    % -- find K nearest neighbours in dictwords of feats --
    if (nargin < 5) || isempty(kdtree)
        distsq = vl_alldist(dictwords, feats);
        
        % distances is MxN matrix where M is num of codewords
        % and N is number of descriptors in feats
        [distsq, ix] = sort(distsq, 1);
        
        % ix is a KxN matrix containing
        % the indices of the K nearest neighbours of each image descriptor
        ix = ix(1:K, :);
        distsq = distsq(1:K, :);        
        
    else
        [ix, distsq] = vl_kdtreequery(kdtree, single(dictwords), ...
                                            single(feats), ...
                                            'MaxComparisons', ...
                                            maxComparisons, ...
                                            'NumNeighbors', K);
        ix = double(ix);
    end
    
    kerMultiplier = 1/(sqrt(2*pi)*sigma);
    kerIntMultiplier = -0.5/(sigma^2);

    % kerDists is a KxN matrix containing the kernel distances of the K
    % nearest neighbours of each image descriptor
    kerDists = kerMultiplier * exp(kerIntMultiplier * distsq);
    
    if isequal(kcbType, 'unc')
        % L1-normalise to ensure soft assignemnts sum to 1
        kerDists = bsxfun(@times, kerDists, 1 ./ max(sum(kerDists, 1), eps));
    end
    
    if outputFullMat        
        % encoding is MxN sparse matrix of results
        % (where M is number of dictionary words, and N is number of feats)
        feat_idx = repmat(1:num_feats, K, 1);
                        
        % construct sparse encoding matrix
        encoding = sparse(ix(:), feat_idx(:), double(kerDists(:)), num_words, num_feats);
    else
        encoding = vl_binsum(zeros(num_words, 1), double(kerDists(:)), ix(:));        
    end    
    
end

