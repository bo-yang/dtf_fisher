classdef GmmCodebkGen < handle & featpipem.codebkgen.GenericCodebkGen
    %GMMCODEBKGEN Generate codebook of visual words using kmeans

    properties
        gauss_count % number of gaussians
        descount_limit % limit on # features to use for clustering
        trainimage_limit % limit on # images to use for clustering
    end
    
    methods
        function obj = GmmCodebkGen(featextr, gauss_count)
            obj.featextr = featextr;
            obj.gauss_count = gauss_count;
            obj.descount_limit = 10e5;
            obj.trainimage_limit = -1;
        end
        codebook = train(obj, imlist, varargin)
    end
    
end

