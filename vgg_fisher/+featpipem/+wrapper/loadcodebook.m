function codebook = loadcodebook(codebkgen, prms)
%LOADCODEBOOK Summary of this function goes here
%   Detailed explanation goes here

if exist(prms.codebook,'file')
    load(prms.codebook);
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
    codebook = codebkgen.train(imfiles);
    
    save(prms.codebook,'codebook','codebkgen');
end

end

