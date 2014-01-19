classdef VQEncoder < handle & featpipem.encoding.GenericEncoder
    %VQENCODER Bag-of-word histogram computation using the VQ method (hard assignment)
    
    properties
        norm_type % 'l1' or 'l2'
        max_comps % -1 for exact
        codebook_
    end
    
    properties(SetAccess=protected)
        kdtree_
    end
    
    methods
        function obj = VQEncoder(codebook)
            % set default parameter values
            obj.norm_type = 'l1';
            obj.max_comps = -1;
            
            % setup encoder
            obj.codebook_ = codebook;
            obj.kdtree_ = vl_kdtreebuild(obj.codebook_);
        end
        
        function dim = get_input_dim(obj)
            dim = size(obj.codebook_,1);
        end
        
        function dim = get_output_dim(obj)
            dim = size(obj.codebook_,2);
        end
                
        % compute encoding
        code = encode(obj, feats)
        
        % compute assignent matrix
        assign = get_assignments(obj, feats)
    end
    
end

