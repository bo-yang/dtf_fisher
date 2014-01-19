function res = dstest(prms, codebook, featextr, encoder, pooler, classifier)
%TEST Summary of this function goes here
%   Detailed explanation goes here

% --------------------------------
% Prepare output filenames
% --------------------------------

trainSetStr = [];
for si = 1:length(prms.splits.train)
    trainSetStr = [trainSetStr prms.splits.train{si}]; %#ok<AGROW>
end

testSetStr = [];
for si = 1:length(prms.splits.test)
    testSetStr = [testSetStr prms.splits.test{si}]; %#ok<AGROW>
end

if ~isfield(prms.experiment, 'classif_tag')
    prms.experiment.classif_tag = '';
end

kChunkIndexFile = fullfile(prms.paths.codes, sprintf('%s_chunkindex.mat', prms.experiment.codes_suffix));
kKernelFile = fullfile(prms.paths.compdata, sprintf('%s_%s_K.mat', prms.experiment.name, trainSetStr));
kClassifierFile = fullfile(prms.paths.compdata, sprintf('%s_%s_classifier%s.mat', prms.experiment.name, trainSetStr, prms.experiment.classif_tag));
kResultsFile = fullfile(prms.paths.results, sprintf('%s_%s_results%s.mat', prms.experiment.name, testSetStr, prms.experiment.classif_tag));

% --------------------------------
% Compute Chunks (for all splits)
% --------------------------------
if exist(kChunkIndexFile,'file')
    load(kChunkIndexFile)
else
    chunk_files = featpipem.chunkio.compChunksIMDB(prms, featextr, pooler);
    % save chunk_files to file
    save(kChunkIndexFile, 'chunk_files');
end


train_chunks = cell(1,length(prms.splits.train));
for si = 1:length(prms.splits.train)
    train_chunks{si} = chunk_files(prms.splits.train{si});
end

% --------------------------------
% Compute Kernel (if using a dual classifier)
% --------------------------------
if isa(classifier, 'featpipem.classification.svm.LibSvmDual')
    if exist(kKernelFile,'file')
        load(kKernelFile);
    else
        K = featpipem.chunkio.compKernel(train_chunks);
        % save kernel matrix to file
        save(kKernelFile, 'K');
    end
end

% --------------------------------
% Train Classifier
% --------------------------------
if isa(classifier, 'featpipem.classification.svm.LibSvmDual')
    % ...........................
    % training for svm in dual
    % ...........................
    if exist(kClassifierFile,'file')
        load(kClassifierFile);
        classifier.set_model(model); %#ok<NODEF>
    else
        labels_train = featpipem.utility.getImdbGT(prms.imdb, prms.splits.train, 'concatOutput', true);
        classifier.train(K, labels_train, train_chunks);
        model = classifier.get_model(); %#ok<NASGU>
        save(kClassifierFile,'model');
    end
else
    % ...........................
    % training for svm in primal
    % ...........................
    if exist(kClassifierFile,'file')
        load(kClassifierFile);
        classifier.set_model(model); %#ok<NODEF>
    else
        labels_train = featpipem.utility.getImdbGT(prms.imdb, prms.splits.train, 'concatOutput', true);
        trainvecs = featpipem.chunkio.loadChunksIntoMat(train_chunks);
        classifier.train(trainvecs, labels_train);
        model = classifier.get_model(); %#ok<NASGU>
        save(kClassifierFile,'model');
        clear trainvecs;
    end
end

% --------------------------------
% Test Classifier
% --------------------------------
scoremat = cell(1,length(prms.splits.test));
res = cell(1,length(prms.splits.test));
% apply classifier to all testsets in prms.splits.test
rankingTime = 0;
for si = 1:length(prms.splits.test)
    tRankingTime = tic;
    [scoremat{si}, scoremat{si}] = featpipem.chunkio.testChunks(chunk_files(prms.splits.test{si}), classifier);
    rankingTime = rankingTime + toc(tRankingTime);
    switch prms.experiment.evalusing
        case 'precrec'
            res{si} = featpipem.eval.evalPrecRec(prms.imdb, scoremat{si}, prms.splits.test{si}, prms.experiment.dataset);
        case 'accuracy'
            res{si} = featpipem.eval.evalAccuracy(prms.imdb, scoremat{si}, prms.splits.test{si});
        otherwise
            error('Unknown evaluation method %s', prms.experiment.evalusing);
    end
end
fprintf('Ranking time was: %f seconds\n',rankingTime);

% package results
results.res = res;
results.scoremat = scoremat; %#ok<STRNU>
parameters.prms = prms;
parameters.codebook = codebook;
parameters.featextr = featextr;
parameters.encoder = encoder;
parameters.pooler = pooler;
parameters.classifier = classifier; %#ok<STRNU>
% save results to file
save(kResultsFile, 'results', 'parameters','-v7.3');

end

