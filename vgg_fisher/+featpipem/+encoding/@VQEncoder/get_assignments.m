function assign = get_assignments(obj, feats)
%GET_ASSIGNMENTS Get hard assignment of features
% Returns MxN sparse matrix of assignments, where M is number of words, and
% N is number of faetures

    % Apply encoding ------------------------------------------------------
    
    if obj.max_comps ~= -1
        % using ann...
        codeids = vl_kdtreequery(obj.kdtree_, obj.codebook_, feats, ...
            'MaxComparisons', obj.max_comps);
    else
        % using exact assignment...
        [~, codeids] = min(vl_alldist(obj.codebook_, feats), [], 1);
    end
    
    num_words = size(obj.codebook_, 2);
    num_feats = size(feats, 2);
    
    assign = sparse(double(codeids), 1:num_feats, 1, num_words, num_feats);
    
end

