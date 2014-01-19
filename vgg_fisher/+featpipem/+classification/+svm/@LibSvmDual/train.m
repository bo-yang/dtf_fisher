function train(obj, input, labels, chunk_files)
%TRAIN Testing function for LIBSVM (using dual formulation)
%   Refer to GenericSVM for interface definition
%   'input' in this case is the precomputed kernel matrix
%   the optional input 'train_codes_or_chunk_files' if specified,
%   cause the computed model to be stored as a w vector as in the primal
%   case, else the SVM model from LIBSVM is stored internally directly

    % process and check input arguments -----

    % if the input is just a cell array of strings, nest it in a second
    % level (as this is the form used when accepting multiple input test
    % sets)
    if ischar(chunk_files{1})
        chunk_files = {chunk_files};
    end
    
    K = input;

    if (size(K,1) ~= size(K,2))
        error('Precomputed kernel matrix is not square');
    end
    
    % prepare temporary output model storage variables
    libsvm = cell(1,length(labels));

    % train models for each class in turn
    %for ci = 1:length(labels)
    parfor ci = 1:length(labels)
        fprintf('Training model for class %d of %d...\n', ci, length(labels));
        labels_cls = -ones(1,size(K,1));
        labels_cls(labels{ci}) = 1;
        
        libsvm{ci} = svmtrain(labels_cls', ...
            [(1:size(K,1))', K+eye(size(K,1))*realmin], ...
            sprintf(' -t 4 -c %f', obj.c)); %#ok<PFBNS>
    end
    
    % get w_full vector from the models if necessary
    obj.modelIsWVect = (nargin > 3);
    if obj.modelIsWVect
        % determine output vector size and preallocate model
        ch = load(chunk_files{1}{1});
        featdim = size(ch.chunk, 1);
        clear ch;
        w_pf = zeros(featdim, length(labels));
        b_pf = zeros(1, length(labels));
        % load in SVs
        %for ci = 1:length(labels)
        parfor ci = 1:length(labels)
            fprintf('Calculating w vector for class %d of %d...\n', ci, length(labels));
            w_pf(:,ci) = single(featpipem.chunkio.compWFromKT(chunk_files, libsvm{ci}.SVs, libsvm{ci}.sv_coef));
            b_pf(ci) = -single(libsvm{ci}.rho);
            % invert labels if necessary
            if labels{ci}(1) ~= 1
                w_pf(:,ci) = w_pf(:,ci)*-1;
                b_pf(ci) = b_pf(ci)*-1;
            end
        end
        % transfer temporary variables across to model store
        obj.model.w = w_pf;
        obj.model.b = b_pf;
    else
        obj.model = libsvm;
    end
end

