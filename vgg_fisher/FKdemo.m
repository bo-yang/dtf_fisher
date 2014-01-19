bCrossValSVM = true;
voc_size = 256; % vocabulary size
desc_dim = 80; % descriptor dimensionality after PCA
DataDir = '/data/ken/enceval_validation_2/';

%% initialize experiment parameters
prms.experiment.name = 'FKtest_comb'; % experiment name - prefixed to all output files other than codes
prms.experiment.codes_suffix = 'FKtest_comb'; % string prefixed to codefiles (to allow sharing of codes between multiple experiments)
prms.experiment.classif_tag = ''; % additional string added at end of classifier and results files (useful for runs with different classifier parameters)
prms.imdb = load(fullfile(DataDir,'imdb/imdb-VOC2007.mat')); % IMDB file
prms.codebook = fullfile(DataDir, sprintf('data/codebooks/fkdemo_gmm_%d_%d.mat', voc_size, desc_dim)); % desired location of codebook
prms.dimred = fullfile(DataDir, sprintf('data/dimred/fkdemo_pca_%d.mat', desc_dim)); % desired location of low-dim projection matrix
prms.experiment.dataset = 'VOC2007'; % dataset name - currently only VOC2007 supported
prms.experiment.evalusing = 'precrec'; % evaluation method - currently only precision recall supported

prms.paths.dataset = '/data/common/datasets/';
prms.paths.codes = fullfile(DataDir,'data/codes/'); % path where codefiles should be stored
prms.paths.compdata = fullfile(DataDir,'data/compdata/'); % path where all other compdata (kernel matrices, SVM models etc.) should be stored
prms.paths.results = fullfile(DataDir,'data/results/'); % path where results should be stored

prms.chunkio.chunk_size = 100; % number of encodings to store in single chunk
prms.chunkio.num_workers = max(matlabpool('size'), 1); % number of workers to use when generating chunks

% initialize split parameters
prms.splits.train = {'train', 'val'}; % cell array of IMDB splits to use when training
prms.splits.test = {'test'}; % cell array of IMDB splits to use when testing

% initialize experiment classes
featextr = featpipem.features.PhowExtractor();
featextr.step = 3;

%% train/load dimensionality reduction
if desc_dim ~= 128
    dimred = featpipem.dim_red.PCADimRed(featextr, desc_dim);
    featextr.low_proj = featpipem.wrapper.loaddimred(dimred, prms);
else
    % no dimensionality reduction
    featextr.low_proj = [];
end

%% train/load codebook
codebkgen = featpipem.codebkgen.GMMCodebkGen(featextr, voc_size);
codebkgen.GMM_init = 'rand';

codebook = featpipem.wrapper.loadcodebook(codebkgen, prms);

%% initialize encoder + pooler
encoder = featpipem.encoding.FKEncoder(codebook);
encoder.pnorm = single(0.0);
encoder.alpha = single(1.0);
encoder.grad_weights = false;
encoder.grad_means = true;
encoder.grad_variances = true;

pooler = featpipem.pooling.SPMPooler(encoder);
pooler.subbin_norm_type = 'l2';
pooler.norm_type = 'none';
pooler.pool_type = 'sum';
pooler.kermap = 'hellinger';
pooler.post_norm_type = 'l2';

%% classification
classifier = featpipem.classification.svm.LibSvmDual();

if bCrossValSVM
    prms.splits.train = {'train'};
    prms.splits.test = {'val'};
    c = [1.6 1.8 2 2.2 2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4 4.2 4.4 4.6 4.8 5 5.2 5.4 5.6 5.8 6 6.2 6.4 6.6 6.8 7 7.2 7.4 7.6 7.8];
    maxmap = 0;
    maxci = 1;
    for ci = 1:length(c)
        prms.experiment.classif_tag = sprintf('c%g', c(ci));
        classifier.c = c(ci);
        AP = featpipem.wrapper.dstest(prms, codebook, featextr, encoder, pooler, classifier);
    
        map = mean(AP{1});
        if map > maxmap
            maxci = ci;
            maxmap = map;
        end
    end
    
    prms.splits.train = {'train', 'val'};
    prms.splits.test = {'test'};
    prms.experiment.classif_tag = sprintf('TESTc%f', c(maxci));
    classifier.c = c(maxci);
    AP = featpipem.wrapper.dstest(prms, codebook, featextr, encoder, pooler, classifier);
else
    classifier.c = 1.6;
    AP = featpipem.wrapper.dstest(prms, codebook, featextr, encoder, pooler, classifier);
end
