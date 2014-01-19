classdef KCBEncoder < handle & featpipem.encoding.GenericEncoder
    %KCBENCODER Bag-of-word histogram computation using the KCB method
    
    properties
        norm_type    % 'l1' or 'l2'
        max_comps    % -1 for exact
        num_nn       % number of nearest neighbour bases to assign to
        sigma        % sigma of Gaussian bases
        kcb_type     % 'unc', 'kcb', 'pla' - KCB method, see paper
        codebook_
    end
    
    properties(SetAccess=protected)
        kdtree_
    end
    
    methods
        function obj = KCBEncoder(codebook)
            % set default parameter values
            obj.norm_type = 'l2';
            obj.max_comps = -1;
            obj.num_nn = 5;
            obj.sigma = 1e-4;
            obj.kcb_type = 'unc';
            
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

