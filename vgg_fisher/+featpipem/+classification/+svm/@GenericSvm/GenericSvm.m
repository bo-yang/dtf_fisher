classdef GenericSvm < handle
    %GENERICSVM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % storage for trained svm model
        model
    end
    
    methods
        function set_model(obj, model)
            obj.model = model;
        end
        function model = get_model(obj)
            model = obj.model;
        end
    end
    
    methods(Abstract)
        % Training Function for SVM ---------------------------------------
        % input - a matrix of column features to use for training
        % labels - a cell array of length K classes, with each cell giving
        %   the indexes of input features containing that class
        train(obj, input, labels)
        % Testing Function for SVM ----------------------------------------
        % input - a matrix of column features to test
        % est_label - the index of the class with the highest score
        % scoremat - the full SVM output matrix of size KxN where K is the
        %   number of classes and N is the number of input features
        [est_label, scoremat] = test(obj, input)
    end
    
end

