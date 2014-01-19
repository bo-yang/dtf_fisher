function mat = loadChunksIntoMat(chunk_files)
%LOADCHUNKSINTOMAT Summary of this function goes here
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
        fprintf('  Loading in features from chunk %d of %d...\n', ci, length(chunk_files{si}));
        
        ch = load(chunk_files{si}{ci});
        % apply index offset if required
        ch.index = ch.index + idxoffset;
        % store maxidx for current set (to calculate offset for next set)
        if ch.index(end) > maxidx, maxidx = ch.index(end); end
        
        % if this is first chunkfile, preallocate output matrix
        if (si == 1) && (ci == 1)
            featdim = size(ch.chunk, 1);
            featcount = size(ch.chunk, 2);
            chunkfilecount = 0;
            for i = 1:length(chunk_files)
                chunkfilecount = chunkfilecount + length(chunk_files{i});
            end
            mat = cast(zeros(featdim, featcount*chunkfilecount), class(ch.chunk));
        end
        
        % now copy over the chunk
        mat(:, ch.index) = ch.chunk;
    end
end

% finally, downsize the output matrix if required
if (maxidx > size(mat,2))
    mat = mat(:,1:maxidx);
end

end

