classdef FKEncoder < handle & featpipem.encoding.GenericEncoder
    %FKENCODER Bag-of-word histogram computation using the FK method
    
    properties
        grad_weights
        grad_means
        grad_variances
        alpha
        pnorm
    end
    
    properties(SetAccess=protected)
        codebook_
    end
    
    properties(SetAccess=protected, GetAccess=protected)
        fc_
        fisher_params_
    end
    
    methods
        function obj = FKEncoder(codebook)
            % set default parameter values
            
            % "soft" BOW
            obj.grad_weights = false;
            
            % 1st order
            obj.grad_means = true;
            
            % 2nd order
            obj.grad_variances = true;
            
            % power normalization (set to 1 to disable)
            obj.alpha = single(0.5);
            
            % norm regularisation (set to 0 to disable)
            obj.pnorm = single(2);
            
            % setup encoder
            obj.codebook_ = codebook;
        end
        function dim = get_input_dim(obj)
            dim = obj.codebook_.n_dim;
        end
        function dim = get_output_dim(obj)
            dim = 0;
            
            if obj.grad_weights
                dim = dim + obj.codebook_.n_gauss;
            end
            
            if obj.grad_means
                dim = dim + obj.codebook_.n_gauss*obj.codebook_.n_dim;
            end
            
            if obj.grad_variances
                dim = dim + obj.codebook_.n_gauss*obj.codebook_.n_dim;
            end
        end
        code = encode(obj, feats)
    end
    
end

