classdef LLCEncoder < handle & featpipem.encoding.GenericEncoder
    %LLCENCODER Bag-of-word histogram computation using the LLC method
    
    properties
        norm_type    % 'l1' or 'l2'
        max_comps    % -1 for exact
        num_nn       % number of nearest neighbour bases to assign to
        beta         % LLC regularization parameter
    end
    
    properties(SetAccess=protected)
        codebook_
        kdtree_
    end
    
    methods
        function obj = LLCEncoder(codebook)
            % set default parameter values
            obj.norm_type = 'l2';
            obj.max_comps = -1;
            obj.num_nn = 5;
            obj.beta = 1e-4;
            
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
        code = encode(obj, feats)
    end
    
end

