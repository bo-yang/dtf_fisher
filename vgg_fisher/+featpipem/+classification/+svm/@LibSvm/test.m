function [est_label, scoremat] = test(obj, input)
%TEST Training function for LIBSVM
%   Refer to GenericSVM for interface definition

    % ensure a model has been trained
    if isempty(obj.model)
        error('A SVM model has yet to be trained');
    end
    
    % ensure input is of correct form
    if ~isa(input,'double')
        input = sparse(double(input));
    end
    
    % prepare output matrix
    scoremat = zeros(length(obj.model.libsvm), size(input,2));
    
    % test models for each class in turn
    for ci = 1:length(obj.model.libsvm)
        [scorevec scorevec scorevec] = ...
            svmpredict(zeros(size(input,2),1), input', ...
            obj.model.libsvm{ci}); %#ok<PFBNS,ASGLU>
        if obj.model.libsvm_flipscore(ci)
            scorevec = -scorevec;
        end
        scoremat(ci,:) = scorevec';
    end
    
    [est_label est_label] = max(scoremat, [], 1); %#ok<ASGLU>
end

