function w = compWFromKT(chunk_files, svs, sv_coefs)
%COMPWFROMKT Summary of this function goes here
%   Detailed explanation goes here

% process and check input arguments -----

if ~iscell(chunk_files)
    error(['chunk_files must be either a cell of strings, or cell ' ...
        'of cell of strings (to support multiple input sets)']);
end
% if the input is just a cell array of strings, nest it in a second
% level (as this is the form used when accepting multiple input test
% sets)
if ischar(chunk_files{1})
    chunk_files = {chunk_files};
end

idxoffset = 0;
maxidx = 0;
% iterate over chunkfiles
for si = 1:length(chunk_files)
    fprintf('Processing set %d of %d...\n', si, length(chunk_files));
    idxoffset = idxoffset+maxidx;
    maxidx = 0;
    for ci = 1:length(chunk_files{si})
        fprintf('  Adding support vectors from chunk %d of %d...\n', ci, length(chunk_files{si}));

        ch = load(chunk_files{si}{ci});
        % apply index offset if required
        ch.index = ch.index + idxoffset;
        % store maxidx for current set (to calculate offset for next set)
        if ch.index(end) > maxidx, maxidx = ch.index(end); end
        % if this is first chunkfile, preallocate output matrix
        if (si == 1) && (ci == 1)
            featdim = size(ch.chunk, 1);
            w = zeros(featdim, 1);
        end

        % for the current datafile, get all support vectors which are
        % contained within it
        [where loc] = ismember(svs,ch.index);
        loc(loc == 0) = [];
        % if there are no support vectors in the current datafile, skip
        % to next one
        if isempty(loc)
            continue;
        end
        % get the w so far by multiplying each of the support vectors by
        % the support vector coefficients, then summing. nFeatures is the
        % size of a single encoding, so we are doing:
        % [featsz x #svcoefs] *. [featsz x #svcoefs]
        W = ch.chunk(:,loc).*repmat(sv_coefs(where)',[featdim 1]);
        w = w+sum(W,2);
    end
end

end

