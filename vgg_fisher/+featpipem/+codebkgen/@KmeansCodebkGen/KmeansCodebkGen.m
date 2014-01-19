classdef KmeansCodebkGen < handle & featpipem.codebkgen.GenericCodebkGen
    %KMEANSCODEBKGEN Generate codebook of visual words using kmeans
    
    properties
        cluster_count % number of visual words in codebook
        descount_limit % limit on # features to use for clustering
        trainimage_limit % limit on # images to use for clustering
        maxcomps % maximum number of comparisons when using ANN (-1 = exact)
    end
    
    methods
        function obj = KmeansCodebkGen(featextr, cluster_count)
            obj.featextr = featextr;
            obj.cluster_count = cluster_count;
            obj.descount_limit = 1e6;
            obj.trainimage_limit = -1;
            obj.maxcomps = ceil(cluster_count/4);
        end
        codebook = train(obj, imlist, varargin)
    end
    
end
