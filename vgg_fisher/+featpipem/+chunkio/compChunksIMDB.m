function chunk_files = compChunksIMDB(prms, featextr, encpooler)
%COMPCHUNKSIMDB Compute feature chunks, save to disk and return filenames
%   Given a set of test parameters 'prms', computes the features for all
%   image sets in a given imdb and saves to chunk files
%
%   chunk_files - output filenames. Is an instance of containers.map, with
%                 one keyed value for each of 'train', 'test' and (if it
%                 exists) 'val', each containing a cell array of chunk
%                 filenames

% default parameters ------------------------------------------------------
% imdb
% paths.dataset
% paths.codes
% paths.compdata
% paths.results
% experiment.name
% experiment.dataset
% experiment.codes_suffix
% chunkio.chunk_size
% chunkio.num_workers
% 

% output format of chunk filenames
% placeholders: codes_suffix targetset chunkstartidx
CHUNK_FNAME = '%s_%s_chunk%d.mat';

% initialize output map
chunk_files = containers.Map();

% iterate over sets in IMDB
for targetsets = {'train', 'val', 'test'}
    
    targetset = targetsets{1};
    switch targetset
        case 'train'
            ids = find(prms.imdb.images.set == prms.imdb.sets.TRAIN);
        case 'val'
            if isfield(prms.imdb.sets, 'VAL')
                ids = find(prms.imdb.images.set == prms.imdb.sets.VAL);
            else
                disp('Skipping validation set (deson''t exist in IMDB)');
                continue;
            end
        case 'test'
            ids = find(prms.imdb.images.set == prms.imdb.sets.TEST);
    end
    
    % calculate chunk start indexes
    chunk_starts_ = 1:prms.chunkio.chunk_size:length(ids);
    
    % allocate chunks to workers
    chunk_starts = cell(prms.chunkio.num_workers);
    for w = 1:prms.chunkio.num_workers
        chunk_starts{w} = chunk_starts_(w:prms.chunkio.num_workers:end);
    end
    
    % initialize storage for chunk filenames in parfor
    chunk_files_by_worker = cell(prms.chunkio.num_workers, 1);
    
    % compute chunks
%     for w = 1:prms.chunkio.num_workers
    parfor w = 1:prms.chunkio.num_workers
        % iterate through chunks assigned to current worker
        for c = 1:length(chunk_starts{w})
            fprintf('Processing chunk starting at %d (worker: %d round: %d)...\n',...
                chunk_starts{w}(c), w, c);
            path = fullfile(prms.paths.codes, sprintf(CHUNK_FNAME, prms.experiment.codes_suffix, targetset, chunk_starts{w}(c)));
            if exist(path,'file')
                fprintf('Chunkfile exists. Skipping...\n');
                continue;
            end
            this_chunk_size = min(prms.chunkio.chunk_size, length(ids)-chunk_starts{w}(c)+1);
            chunk = zeros(encpooler.get_output_dim, this_chunk_size, 'single');
            % iterate through images in current chunk
            for i = chunk_starts{w}(c):chunk_starts{w}(c)+this_chunk_size-1
                fprintf('  computing features for image: %s....\n', prms.imdb.images.name{ids(i)}); 
                im = imread(fullfile(prms.paths.dataset, prms.imdb.images.name{ids(i)}));
                im = featpipem.utility.standardizeImage(im);
                [feats, frames] = featextr.compute(im);
                cur_code = encpooler.compute(size(im), feats, frames);
                
                if any(isnan(cur_code))
                    error('NaNs in the code of chunk %d (worker: %d round: %d)!', i, w, c);
                end
                
                chunk(:,i-chunk_starts{w}(c)+1) = cur_code;
            end
            index = chunk_starts{w}(c):chunk_starts{w}(c)+this_chunk_size-1;
            save_chunk_(path, chunk, index);
            % append filename to output
            chunk_files_by_worker{w}{end+1} = path;
        end
    end
    
    % copy chunk filenames over to output map
    chunk_files(targetset) = {};
    for w = 1:prms.chunkio.num_workers
        chunk_files(targetset) = [chunk_files(targetset) chunk_files_by_worker{w}];
    end
    
end

fprintf('Features computed!\n');

end

function save_chunk_(filename, chunk, index) %#ok<INUSD>
save(filename,'chunk','index','-v7.3');
end


