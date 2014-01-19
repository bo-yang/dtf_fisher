classdef SPMPooler < handle & featpipem.pooling.GenericPooler
    %SPMPOOLER Pooling using the spatial pyramid match kernel
    
    properties
        subbin_norm_type    % 'l1' or 'l2' (or other value = none)
        norm_type    % 'l1' or 'l2' (or other value = none)
        post_norm_type    % 'l1' or 'l2' (or other value = none)
        pool_type    % 'sum' or 'max'
        quad_divs
        horiz_divs
        kermap  % 'homker', 'hellinger' (or other value = none [default])
    end
    
    properties(SetAccess=protected)
        encoder_     % implementation of featpipem.encoding.GenericEncoder
    end
    
    methods
        function obj = SPMPooler(encoder)
            % set default parameter values
            obj.subbin_norm_type = 'l2';
            obj.norm_type = 'l2';
            obj.post_norm_type = 'none';
            obj.pool_type = 'sum';
            obj.quad_divs = 2;
            obj.horiz_divs = 3;
            obj.kermap = 'none';
            
            % setup encoder
            obj.encoder_ = encoder;
        end
        function dim = get_output_dim(obj)
            % check bin levels
            if mod(log2(obj.quad_divs),1)
                error('quad_divs must be a power of 2');
            end
            bin_quads_count = obj.quad_divs*obj.quad_divs;
            bin_quad_levels = 1;
            bin_div_tmp = obj.quad_divs;
            while bin_div_tmp ~= 2
                bin_div_tmp = bin_div_tmp/2;
                bin_quad_levels = bin_quad_levels + 1;
                bin_quads_count = bin_quads_count + bin_div_tmp*bin_div_tmp;
            end
            clear bin_div_tmp;

            bin_count = bin_quads_count + obj.horiz_divs + 1;
            dim = bin_count*obj.encoder_.get_output_dim();
            % account for expansion in dimensionality when using kernel map
            if strcmp(obj.kermap,'homker')
                dim = dim*3;
            end
        end
        pcode = compute(obj, imsize, feats, frames)
    end
    
end

