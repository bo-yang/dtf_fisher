function [est_label, scoremat] = testChunks(chunk_files, svm)
%TESTCHUNKS Summary of this function goes here
%   Detailed explanation goes here

for c1 = 1:length(chunk_files)        
    % apply classifier
    fprintf('Applying classifier to chunk %d of %d...\n', c1, length(chunk_files));
    ch = load(chunk_files{c1});

    [est_label_ch scoremat_ch] = svm.test(ch.chunk);
    assert(length(est_label_ch) == length(ch.index));

    % preallocate output matrices on first iteration
    if c1 == 1
        est_label = zeros(1,size(est_label_ch,2)*length(chunk_files));
        indices = zeros(1,size(est_label_ch,2)*length(chunk_files));
        scoremat = zeros(size(scoremat_ch,1),size(est_label_ch,2)*length(chunk_files));
        startidx = 1;
    end

    % copy output to output arrays
    maxidx = startidx + length(est_label_ch) - 1;
    est_label(:,startidx:maxidx) = est_label_ch;
    indices(:,startidx:maxidx) = ch.index;
    scoremat(:,startidx:maxidx) = scoremat_ch;
    startidx = maxidx + 1;
end

% downsize output arrays if required
if maxidx < size(est_label,2)
    est_label = est_label(:,1:maxidx);
    indices = indices(:,1:maxidx);
    scoremat = scoremat(:,1:maxidx);
end

% sort in order of indices, and ensure all indices are present
[indices sortidx] = sort(indices);
est_label = est_label(:,sortidx);
scoremat = scoremat(:,sortidx);

index_test = 1:length(indices);
if ~isequal(indices, index_test)
    error('Indices extracted from chunk files are non-contiguous');
end

end

