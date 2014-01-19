function code = encode(obj, feats)
%ENCODE Encode features using the VQ method (hard assignment)

    % Apply encoding ------------------------------------------------------
    
    if obj.max_comps ~= -1
        % using ann...
        codeids = vl_kdtreequery(obj.kdtree_, obj.codebook_, feats, ...
            'MaxComparisons', obj.max_comps);
    else
        % using exact assignment...
        [~, codeids] = min(vl_alldist(obj.codebook_, feats), [], 1);
    end
    
    code = vl_binsum(zeros(size(obj.codebook_, 2), 1), 1, double(codeids));
    code = single(code);
    
    % Normalize -----------------------------------------------------------
    
    if strcmp(obj.norm_type, 'l1')
        code = code / norm(code,1);
    end
    if strcmp(obj.norm_type, 'l2')
        code = code / norm(code,2);
    end

end

