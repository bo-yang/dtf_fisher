classdef GenericCodebkGen < handle
    %GENERICCODEBKGEN Generic interface for training codebooks
    
    properties
        featextr
    end
    
    methods(Abstract)
        codebook = train(obj, imlist)
    end
    
end

