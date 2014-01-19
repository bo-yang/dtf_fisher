classdef PCADimRed < handle & featpipem.dim_red.GenericDimRed
    %PCADimRed Learns descriptor dimensionality reduction using PCA
    
    properties
        dim % target dimensionality
        descount_limit % limit on # features to use
        trainimage_limit % limit on # images to use        
    end
    
    methods
        function obj = PCADimRed(featextr, dim)
            
            obj.featextr = featextr;
            obj.dim = dim;
            obj.descount_limit = 1e6;
            obj.trainimage_limit = -1;            
        end
        
        low_proj = train(obj, imlist, varargin)
    end
    
end

