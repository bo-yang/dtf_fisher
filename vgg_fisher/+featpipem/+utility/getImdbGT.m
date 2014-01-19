function [gt set_sizes] = getImdbGT(imdb, sets, varargin)
%GETIMDBGT Get ground truth from IMDB file
%
%   Returns ground truth for one or more named sets from an IMDB
%
%   INPUTS:
%   imdb - the IMDB input file to use
%   sets - a cell array of sets e.g. {'train', 'val', 'test'}
%   OPTIONAL VARARGIN INPUTS:
%   outputSignedLabels - (optional) toggle format of output gt (see below)
%   concatOutput - (optional) concatenate labels into a single output
%                  vector (see below)
%
%   OUTPUTS:
%   gt - a cell array of size CxS where C is the number of classes in the
%        IMDB and S is the number of sets specified in 'sets'. Each cell is
%        either:
%        IF outputSignedLabels is FALSE (default):
%           a vector of 1-indexed indices for each set, specifying positive
%           samples
%        IF outputSignedLabels is TRUE:
%           a vector of length (# of samples in set) with -1 and +1 values
%           indicating negative and positive samples respectively
%   set_sizes - a vector of length (# sets) specifying the number of
%               samples in each set
%
%   IF concatOutput is TRUE, then instead of a CxS cell array, gt is a Cx1
%   cell array with the output from all sets concatenated into a single
%   vector. If outputSignedLabels is FALSE, the indeces of the second set
%   onwards will be incremented by the total number of samples in the sets
%   before it, so indeces are relative to a a concatenated list comprising
%   all sets.
%
%   If no positive samples are found for any class/set combination, an
%   error is raised

% set default arguments
opts.outputSignedLabels = false;
opts.concatOutput = false;
opts = vl_argparse(opts, varargin);

% get number of classes
n_classes = length(imdb.classes.name);

% prepare output arrays
gt = cell(n_classes, length(sets));
set_sizes = zeros(1, length(sets));

for ci = 1:n_classes
    for si = 1:length(sets)
        % search for set in IMDB
        if ~isfield(imdb.sets, upper(sets{si}))
            error(['Set ''' sets{si} ''' could not be found in IMDB']);
        end
        setid = imdb.sets.(upper(sets{si}));

        % get ids of images in testset
        image_ids = imdb.images.id(imdb.images.set == setid);

        % set output vector to INDEX of images in the set which are
        % also contained within the ground truth of the current class
        [gt{ci,si}, gt{ci,si}] = ...
            intersect(image_ids, imdb.classes.imageIds{ci});
        if isempty(gt{ci,si})
            error(['No ground truth found for specified testset: ' sets{si}]);
        end
        
        % fill set_sizes if first iteration
        if ci == 1
            set_sizes(si) = sum(imdb.images.set == setid);
        end
    end
end

% convert output to signed labels if required
if opts.outputSignedLabels
    for ci = 1:n_classes
        for si = 1:length(sets)
            gt_tmp = -1*ones(1,set_sizes(si));
            gt_tmp(gt{ci,si}) = 1;
            gt{ci,si} = gt_tmp;
        end
    end
end

% convert output to concatenated output if required
if opts.concatOutput
    % deal with indexed output case
    if ~opts.outputSignedLabels
        % add offset indexes to the second set onwards
        for si = 2:length(set_sizes)
            for ci = 1:size(gt,1)
                gt{ci,si} = gt{ci,si} + sum(set_sizes(1:(si-1)));
            end
        end
    end
    % concatenate along set dimension
    gt_final = cell(size(gt,1),1);
    for ci = 1:size(gt,1)
        gt_final{ci} = cell2mat(gt(ci,:));
    end
    gt = gt_final;
end

end

