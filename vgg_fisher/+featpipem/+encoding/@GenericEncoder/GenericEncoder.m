classdef GenericEncoder < handle
    %GENERICENCODER Generic interface to bag-of-words encoder
    
    properties
    end
    
    methods(Abstract)
        dim = get_input_dim(obj)
        dim = get_output_dim(obj)
        code = encode(obj, feats)
    end
    
end

