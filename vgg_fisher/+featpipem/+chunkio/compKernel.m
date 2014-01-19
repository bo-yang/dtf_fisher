function K = compKernel(chunk_files)
%COMPKERNEL Compute kernel matrix from feature chunk files
%   chunk_files -  cell array of paths to chunk files to use to compute K,
%                  with assumed complete monotonically increasing set of
%                  indeces

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

% preallocate kernel matrix
ch = load(chunk_files{1}{1});

size_chunk = size(ch.chunk, 2);
size_est = 0;

for i = 1:length(chunk_files)
    size_est = size_est + size_chunk * length(chunk_files{i});
end

K = zeros(size_est);
clear ch;

idxoffseti = 0;
maxidxi = 0;
% keep track of largest index stored in kernel matrix (used for
% preallocation)
maxidx_ker = 0;

% do_norm = false;

% iterate over first chunkfile
for si = 1:length(chunk_files)
    
    idxoffseti = idxoffseti+maxidxi;
        
%     for ci = 1:length(chunk_files{si})
    parfor ci = 1:length(chunk_files{si})
        
        fprintf('Computing datafile for chunk %d of %d (in set %d of %d)\n', ci, length(chunk_files{si}), si, length(chunk_files));
        
        % load chunk
        ch1 = load(chunk_files{si}{ci});
        
%         if do_norm
%             ch1.chunk = bsxfun(@times, ch1.chunk, 1 ./ sqrt(sum(ch1.chunk .^ 2, 1)));
%         end
        
        % apply index offset if required
        idx1{ci} = ch1.index + idxoffseti;
        
        % part of K
        K1{ci} = zeros(size(ch1.chunk, 2), size_est);
        
        % iterate over second chunkfile
        idxoffsetj = 0;
        maxidxj = 0;
        
        for sj = 1:length(chunk_files)
            
            idxoffsetj = idxoffsetj + maxidxj;
            
            idx2end = zeros(length(chunk_files{sj}), 1);
            
            for cj = 1:length(chunk_files{sj})
                
                ch2 = load(chunk_files{sj}{cj});
                
%                 if do_norm
%                     ch2.chunk = bsxfun(@times, ch2.chunk, 1 ./ sqrt(sum(ch2.chunk .^ 2, 1)));
%                 end
                
                % apply index offset if required
                ch2.index = ch2.index + idxoffsetj;
                idx2end(cj) = ch2.index(end);
                
                % do computation of sub-part of kernel matrix                                
                K1{ci}(:, ch2.index) = ch1.chunk'*ch2.chunk;
                
            end
            
            % store maxidxj for current set (to calculate offset for next set)
            maxidxj = max(idx2end);
        end
    end
    
    % copy the data from K1 to K
    for ci = 1:length(chunk_files{si})
        K(idx1{ci}, :) = K1{ci};
    end
    
    % store maxidxi for current set (to calculate offset for next set)
    maxidxi = max(cellfun(@(x) x(end), idx1));
    
    % store absolute max index to aid resizing of kernel matrix at end
    maxidx_ker = max(maxidx_ker, maxidxi);
end

% use maxidx_ker to resize output matrix to be correct size
K = K(1:maxidx_ker, 1:maxidx_ker);

fprintf('Kernel matrix computed\n');

end

