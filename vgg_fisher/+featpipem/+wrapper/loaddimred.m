function low_proj = loaddimred(dimred, prms)
%loaddimred Summary of this function goes here
%   Detailed explanation goes here

if exist(prms.dimred,'file')
    load(prms.dimred);
else
    trainval_files = prms.imdb.images.name(...
                    prms.imdb.images.set == prms.imdb.sets.TRAIN | ...
                    prms.imdb.images.set == prms.imdb.sets.VAL);
                
    num_files = numel(trainval_files);
    
    imfiles = cell(num_files, 1);
    
    for i = 1:num_files
        imfiles{i} = fullfile(prms.paths.dataset, trainval_files{i});
    end
    
    % do training...
    low_proj = dimred.train(imfiles);
    
    save(prms.dimred, 'low_proj');
end

end

