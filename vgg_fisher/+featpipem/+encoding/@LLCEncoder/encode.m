function code = encode(obj, feats)
%ENCODE Encode features using the LLC method

    % Apply encoding ------------------------------------------------------
    
    if obj.max_comps ~= -1
        % using ann...
        code = featpipem.lib.LLCEncode(feats, obj.codebook_, ...
            obj.num_nn, obj.beta, obj.kdtree_, obj.max_comps, false); %#ok<ASGLU>
    else
        % using exact assignment...
        code = featpipem.lib.LLCEncode(feats, obj.codebook_, ...
            obj.num_nn, obj.beta, [], [], false); %#ok<ASGLU>
    end
    
    % Normalize -----------------------------------------------------------
    
    if strcmp(obj.norm_type, 'l1')
        code = code / norm(code,1);
    end
    if strcmp(obj.norm_type, 'l2')
        code = code / norm(code,2);
    end

end

