classdef LTMatrix < handle
    %LTMATRIX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    properties(SetAccess=protected)
        N
        matstore_
    end
    
    methods
        function obj = LTMatrix(N)
            obj.N = N;
            % calculate the length of storage required to store LTMatrix
            ltmatlen = 0;
            for coladd = 1:(N-1)
                ltmatlen = ltmatlen + coladd;
            end
            % initialize matrix
            obj.matstore_ = zeros(ltmatlen,1);
        end
        val = readIndex(obj, i, j)
        val = writeIndex(obj, i, j, val)
        mat = convToFullMatrix(obj)
        
        linidx = getIndex_(obj, i, j)
    end
    
end

