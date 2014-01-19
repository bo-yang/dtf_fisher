function WMat = getWMat(obj)
%GETWMAT Summary of this function goes here
%   Detailed explanation goes here

    % ensure a model has been trained
    if isempty(obj.model)
        error('A SVM model has yet to be trained');
    end
    
    WMat = zeros(size(obj.model.libsvm{1}.SVs,2)+1, length(obj.model.libsvm));
    
    for ci = 1:length(obj.model.libsvm)
        WMat(1:(end-1),ci) = obj.model.libsvm{ci}.SVs'*obj.model.libsvm{ci}.sv_coef;
        WMat(end,ci) = -obj.model.libsvm{ci}.rho;
        if obj.model.libsvm_flipscore{ci}
            WMat(:,ci) = -WMat(:,ci);
        end
    end
end

