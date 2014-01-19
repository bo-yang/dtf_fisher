classdef GenericPooler < handle
    %GENERICPOOLER Generic interface to pooler
    
    properties
    end
    
    methods(Abstract)
        dim = get_output_dim(obj)
        pcode = compute(obj, imsize, feats, frames)
    end
    
end

