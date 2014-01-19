classdef LinearSvm < handle & featpipem.classification.svm.GenericSvm
    %LINEARSVM Summary of this class goes here
    %   Detailed explanation goes here
    
    methods(Abstract)
        WMat = getWMat(obj);
    end
    
end

