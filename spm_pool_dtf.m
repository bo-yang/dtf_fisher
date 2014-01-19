function pcode = spm_pool_dtf(params, codebook, imsize, feats, frames)
%COMPUTE Pool features using the spatial pyramid match kernel

    % check pool type is valid
    if ~strcmp(params.pool_type,'sum') && ~strcmp(params.pool_type,'max')
        error('pool_type must be either ''sum'' or ''max''');
    end

    % check bin levels
    if mod(log2(params.quad_divs),1)
        error('quad_divs must be a power of 2');
    end
    bin_quads_count = params.quad_divs*params.quad_divs;
    bin_quad_levels = 1;
    bin_div_tmp = params.quad_divs;
    while bin_div_tmp ~= 2
        bin_div_tmp = bin_div_tmp/2;
        bin_quad_levels = bin_quad_levels + 1;
        bin_quads_count = bin_quads_count + bin_div_tmp*bin_div_tmp;
    end
    clear bin_div_tmp;
    
    bin_count = bin_quads_count + params.horiz_divs + 1;
    
    pcode = zeros(get_output_dim(params,codebook), bin_count, 'single');
    
    % first compute for finest 'quarter' bins
    h_unit = imsize(2) / params.quad_divs;
    w_unit = imsize(1) / params.quad_divs;
    y_bin = ceil(frames(2,:) / h_unit);
    x_bin = ceil(frames(1,:) / w_unit);
    
    feats_sel_num = zeros(params.quad_divs ^ 2, 1);
    code_idx = 0;
    
    for sx_bin = 1:params.quad_divs
        for sy_bin = 1:params.quad_divs
            
            code_idx = code_idx + 1;
            
            feats_sel = feats(:, (y_bin == sy_bin) & (x_bin == sx_bin));
            feats_sel_num(code_idx) = size(feats_sel, 2);
            
            if ~isempty(feats_sel)
                %TO FIX
                code = fisher_encode(feats_sel,codebook, params);
                if nnz(isnan(code)), error('Code contains NaNs'); end
                pcode(:, code_idx) = code;
            else
                warning('SPMPool:EmptyBin','empty bin!');
            end
        end
    end
    
    % now merge to form subsequent bin levels
    if bin_quad_levels > 1
        prev_level_bin_divs = params.quad_divs;
        prev_level_start_bin = 1;
        for i = 2:bin_quad_levels
            level_bin_divs = prev_level_bin_divs/2;
            level_start_bin = prev_level_start_bin + ...
                prev_level_bin_divs*prev_level_bin_divs;
            
            for sx_bin = 1:level_bin_divs
                for sy_bin = 1:level_bin_divsobj
                    xstart = 2*sx_bin-1;
                    xend = 2*(sx_bin+1)-2;
                    ystart = 2*sy_bin-1;
                    yend = 2*(sy_bin+1)-2;
                    bins_sel = zeros((xend-xstart+1)*(yend-ystart+1),1);
                    bsidx = 1;
                    for xri = xstart:xend
                        for yri = ystart:yend
                            bins_sel(bsidx) = prev_level_start_bin + ...
                                prev_level_bin_divs*(xri-1) + yri - 1;
                            bsidx = bsidx + 1;
                        end
                    end
                    
                    % combine bins
                    if strcmp(params.pool_type,'sum')
                        pcode(:,level_start_bin + ...
                            level_bin_divs*(sx_bin-1) + ...
                            sy_bin - 1) = sum(pcode(:,bins_sel),2);
                    end
                    if strcmp(params.pool_type,'max')
                        pcode(:,level_start_bin + ...
                            level_bin_divs*(sx_bin-1) + ...
                            sy_bin - 1) = max(pcode(:,bins_sel),[],2);
                    end
                end
            end
            
            prev_level_start_bin = level_start_bin;
            prev_level_bin_divs = level_bin_divs;
        end
    end
    
    % now compute the three horizontal bins
    if params.horiz_divs > 0
        h_hunit = imsize(2) / params.horiz_divs;
        h_ybin = ceil(frames(2,:) / h_hunit);
        for sy_bin = 1:params.horiz_divs
            feats_sel = feats(:, (h_ybin == sy_bin));
            if ~isempty(feats_sel)
                % FIX THIS !!!
                code =  fisher_encode(feats_sel,codebook, params);
%                 if nnz(isnan(code)), error('Code contains NaNs'); end
                pcode(:,bin_quads_count+sy_bin) = code;
            else
                warning('SPMPool:EmptyBin','empty bin!');
            end
        end
    end
    
    if true
        
        % now combine the first binset to make up the whole image bin
        if strcmp(params.pool_type,'sum')
            pcode(:, end) = sum(pcode(:, 1 : params.quad_divs ^ 2),2);% * (feats_sel_num / sum(feats_sel_num));
        elseif strcmp(params.pool_type,'max')
            pcode(:,end) = max(pcode(:,1:params.quad_divs^2),[],2);
        end
        
    else    
        % compute whole image encoding from scratch
        % FIX THIS !!!
        pcode(:, end) = fisher_encode(feats,codebook, obj);    
    end
    
%     if nnz(pcode(:,end)) < 1, error('Code is all zeros!'); end
    
    % now normalize all sub-bins
    if strcmp(params.subbin_norm_type, 'l2')
        
        pcode_norm = sqrt(sum(pcode .^ 2, 1));
        pcode_norm = max(pcode_norm, eps);
        pcode = bsxfun(@times, pcode, 1 ./ pcode_norm);
        
    elseif strcmp(params.subbin_norm_type, 'l1')
        
        pcode_norm = sum(pcode, 1);
        pcode_norm = max(pcode_norm, eps);
        pcode = bsxfun(@times, pcode, 1 ./ pcode_norm);
        
    end
    
    % vectorise
    pcode = pcode(:);
    
    % now normalize whole code
    if strcmp(params.norm_type,'l2')
        pcode = pcode/norm(pcode,2);
    elseif strcmp(params.norm_type,'l1')
        pcode = pcode/norm(pcode,1);
    end
    
    % now apply kernel map if specified
    if ~isequal(params.kermap, 'none')        
        % (note: when adding extra kernel maps, note that the getDim function
        % must also be modified to reflect the appropriate increase in code
        % dimensionality)
        if strcmp(params.kermap,'homker')
            % chi-squared approximation
            pcode = vl_homkermap(pcode, 1, 'kchi2');
        elseif strcmp(params.kermap,'hellinger')
            % "generalised" (signed) Hellinger kernel
            pcode = sign(pcode) .* sqrt(abs(pcode));        
        end

        % now post-normalize whole code
        if strcmp(params.post_norm_type,'l2')
            pcode = pcode/norm(pcode,2);
        elseif strcmp(params.post_norm_type,'l1')
            pcode = pcode/norm(pcode,1);
        end
    end
end


 function dim = get_output_dim(obj,codebook)
            dim = 0;
            
            if obj.grad_weights
                dim = dim + codebook.n_gauss;
            end
            
            if obj.grad_means
                dim = dim + codebook.n_gauss*codebook.n_dim;
            end
            
            if obj.grad_variances
                dim = dim + codebook.n_gauss*codebook.n_dim;
            end
end

