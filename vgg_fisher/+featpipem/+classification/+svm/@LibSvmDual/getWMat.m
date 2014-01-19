function WMat = getWMat(obj)
%GETWMAT Summary of this function goes here
%   Detailed explanation goes here

    % ensure a model has been trained
    if isempty(obj.model)
        error('A SVM model has yet to be trained');
    end
    
    if ~obj.modelIsWVect
        error('SVM model must be W Vector to use this function');
    end
    
    WMat = [obj.model.w; ones(1, size(obj.model.w,2)).*obj.model.b];
end

