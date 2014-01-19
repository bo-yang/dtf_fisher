function code = encode(obj, feats)
%ENCODE Encode features using the KCB method

    % Apply encoding ------------------------------------------------------
    
    if obj.max_comps ~= 1
        % using ann...
        code = featpipem.lib.KCBEncode(feats, obj.codebook_, obj.num_nn, ...
            obj.sigma, obj.kdtree_, obj.max_comps, obj.kcb_type, false);
    else
        % using exact assignment...
        code = featpipem.lib.KCBEncode(feats, obj.codebook_, obj.num_nn, ...
            obj.sigma, [], [], obj.kcb_type, false);
    end
    
    % Normalize -----------------------------------------------------------
    
    if strcmp(obj.norm_type, 'l1')
        code = code / norm(code,1);
    end
    if strcmp(obj.norm_type, 'l2')
        code = code / norm(code,2);
    end

end

