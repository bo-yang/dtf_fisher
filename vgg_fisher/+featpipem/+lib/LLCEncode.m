function encoding = LLCEncode(imwords, dictwords, K, beta, kdtree, ...
    maxComparisons, outputFullMat)


    % START validate input parameters -------------------------------------
    default('K', 5);
    default('beta', 1e-4);
    % only used if a kdtree is specified
    default('maxComparisons', size(dictwords,2));
    default('outputFullMat', true);
    % END validate input parameters ---------------------------------------

    % -- find K nearest neighbours in dictwords of imwords --
    if (nargin < 5) || isempty(kdtree)
        distances = vl_alldist2(double(dictwords),double(imwords));
        % distances is MxN matrix where M is num of codewords
        % and N is number of descriptors in imwords
        [ix, ix] = sort(distances); %#ok<ASGLU>
        % ix is a KxN matrix containing
        % the indices of the K nearest neighbours of each image descriptor
        ix(K+1:end,:) = [];
    else
        ix = vl_kdtreequery(kdtree, single(dictwords), ...
                                            single(imwords), ...
                                            'MaxComparisons', ...
                                            maxComparisons, ...
                                            'NumNeighbors', K);
    end
    
    encoding = ...
            featpipem.lib.LLCEncodeHelper(double(dictwords), double(imwords), ...
                                          double(ix), double(beta), outputFullMat);
    

    % Below is the implementation provided by the Matlab sources
    % available by the paper authors at:
    % http://www.ifp.illinois.edu/~jyang29/LLC.htm

    % % find k nearest neighbors
    % XX = sum(X.*X, 2);
    % BB = sum(B.*B, 2);
    % D  = repmat(XX, 1, nbase)-2*X*B'+repmat(BB', nframe, 1);
    % IDX = zeros(nframe, knn);
    % for i = 1:nframe,
    % 	d = D(i,:);
    % 	[dummy, idx] = sort(d, 'ascend');
    % 	IDX(i, :) = idx(1:knn);
    % end
    % 
    % % llc approximation coding
    % II = eye(knn, knn);
    % Coeff = zeros(nframe, nbase);
    % for i=1:nframe
    %    idx = IDX(i,:);
    %    z = B(idx,:) - repmat(X(i,:), knn, 1);   % shift ith pt to origin
    %    C = z*z';                                % local covariance
    %    C = C + II*beta*trace(C);                % regularlization (K>D)
    %    w = C\ones(knn,1);
    %    w = w/sum(w);                            % enforce sum(w)=1
    %    Coeff(i,idx) = w';
    % end
end
