classdef GenericDimRed < handle
    %GenericDimRed Generic interface for dimensionality reduction
    
    properties
        featextr
    end
    
    methods(Abstract)
        low_proj = train(obj, imlist)
    end
    
end

