function assign = get_assignments(obj, feats)
%GET_ASSIGNMENTS Get soft assignment of features
% Returns MxN sparse matrix of assignments, where M is number of words, and
% N is number of faetures

    % Apply encoding ------------------------------------------------------
    
    if obj.max_comps ~= -1
        % using ann...
        assign = featpipem.lib.KCBEncode(feats, obj.codebook_, obj.num_nn, ...
            obj.sigma, obj.kdtree_, obj.max_comps, obj.kcb_type, true);
    else
        % using exact assignment...
        assign = featpipem.lib.KCBEncode(feats, obj.codebook_, obj.num_nn, ...
            obj.sigma, [], [], obj.kcb_type, true);
    end    
    
end

