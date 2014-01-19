function code = encode(obj, feats)
%ENCODE Encode features using the FK method

    % Initialize encoder wrapper if not done already ----------------------

    if isempty(obj.fc_)
        obj.fisher_params_.grad_weights = obj.grad_weights;
        obj.fisher_params_.grad_means = obj.grad_means;
        obj.fisher_params_.grad_variances = obj.grad_variances;
        obj.fisher_params_.alpha = obj.alpha;
        obj.fisher_params_.pnorm = obj.pnorm;
        
        obj.fc_ = FisherEncoder(obj.codebook_, obj.fisher_params_);
    else
        if ((obj.fisher_params_.grad_weights ~= obj.grad_weights) || ...
                (obj.fisher_params_.grad_means ~= obj.grad_means) || ...
                (obj.fisher_params_.grad_variances ~= obj.grad_variances) || ...
                (obj.fisher_params_.alpha ~= obj.alpha) || ...
                (obj.fisher_params_.pnorm ~= obj.pnorm))
            error(['Fisher parameters cannot be ' ...
                'changed between calls to ''encode()''']);
        end
    end

    % Apply encoding ------------------------------------------------------
    
    disp(size(feats));
    tic;
    code = obj.fc_.encode(feats);
    toc;

end

