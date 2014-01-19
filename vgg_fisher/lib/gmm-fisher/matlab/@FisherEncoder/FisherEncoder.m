classdef FisherEncoder < handle
    %FISHERENCODER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Hidden = true, SetAccess = private)
        cpp_handle;
    end
    
    methods
        % constructor
        function this = FisherEncoder(gmm_struct, fisher_params)
            if (nargin < 1)
                error('Must pass in GMM Structure to initialise encoder');
            end
            if (nargin < 2) || isempty(fisher_params)
                this.cpp_handle = mexFisherEncodeHelperSP('init', gmm_struct);
            else
                this.cpp_handle = mexFisherEncodeHelperSP('init', gmm_struct, fisher_params);
            end
        end
        % destructor
        function delete(this)
            mexFisherEncodeHelperSP('clear', this.cpp_handle);
        end
        % encoding method
        function code = encode(this, x, weights)
            if (nargin < 2)
                error('Must pass in matrix of vectors to encode');
            end
            if (nargin < 3) || isempty(weights)
                code = mexFisherEncodeHelperSP('encode', this.cpp_handle, x);
            else
                code = mexFisherEncodeHelperSP('encode', this.cpp_handle, x, weights);
            end
        end
        % get FK dimensionality
        function dim = getdim(this)
            dim = mexFisherEncodeHelperSP('getdim', this.cpp_handle);
        end
    end
    
end

