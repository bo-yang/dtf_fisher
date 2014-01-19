% Script to test DTF features and Fisher Vector

run('vlfeat-0.9.17/toolbox/vl_setup')

path=pwd;
addpath(path);
addpath(fullfile(path,'yael'));
addpath(fullfile(path,'yael/matlab'));
addpath(fullfile(path,'libsvm'));
addpath(fullfile(path,'libsvm/matlab'));
addpath(fullfile(path,'vgg_fisher'));
addpath(fullfile(path,'vgg_fisher/lib/gmm-fisher/matlab'));

num_feats_per_video=1000; % number of DTF features subsampled from each video
encoding='fisher'; % encoding type: 'fisher' - fisher vector; 'bow' - bag-of-words

dtf_dir=fullfile(path,'UCF101_DTF');
train_dataset=fullfile(path,'data',sprintf('UCF101_train_%ddtf.mat',num_feats_per_video));
test_dataset=fullfile(path,'data',sprintf('UCF101_test_%ddtf.mat',num_feats_per_video));
pred_results=fullfile(path,'data','pred_results.mat');
fisher_vectors=fullfile(path,'data','fisher_vectors.mat');

%% Load DTF features for train data
if ~exist(train_dataset,'file')
    % extract DTF data and labels
    train_list_dir=fullfile(path, 'ucfTrainTestlist'); 
    fprintf('Loading DTF trainning features...\n');
    [dtf_train, train_label]=load_dtf_data(dtf_dir,train_list_dir,'train*', encoding, num_feats_per_video);
    save(train_dataset,'dtf_train','train_label','-v7.3');
else
    load(train_dataset);
end

%% Load DTF features for test data
if ~exist(test_dataset,'file')  
    % extract DTF data
    test_list_dir=fullfile(path, 'ucfTrainTestlist'); 
    fprintf('Loading DTF test features...\n');
    [dtf_test,test_label]=load_dtf_data(dtf_dir,test_list_dir,'test*', encoding, num_feats_per_video);
    save(test_dataset,'dtf_test','test_label','-v7.3');
else
    load(test_dataset);
end

switch encoding
      case 'fisher'
        %% Construct Fisher Vector
        fprintf('Constructing Fisher Vectors...\n');
        K=64;
        num_subsamp=K*2000;
        fvt_train=compute_fisher_vector_vgg( dtf_train, K , num_subsamp ,num_feats_per_video );
        fvt_test=compute_fisher_vector_vgg( dtf_test, K , num_subsamp ,num_feats_per_video );
        %fvt_train=compute_fisher_vector_yael( dtf_train, K , num_subsamp ,num_feats_per_video );
        %fvt_test=compute_fisher_vector_yael( dtf_test, K , num_subsamp ,num_feats_per_video );
        save(fisher_vectors,'fvt_train','fvt_test','train_label','test_label','-v7.3');
        
        % SVM classification
        fprintf('Doing SVM classification...\n');
        pred=cell(length(fvt_train),1);
        ac=cell(length(fvt_train),1);
        svm_param='-t 0 -q -c 100'; % Linear SVM
        
        %for i=1:length(pred)
        for i=length(pred):length(pred) %%% TEST ONLY %%%
            model=ovrtrain(train_label',double(fvt_train{i}'),svm_param);
            [pred{i}, ac{i}, ~]=ovrpredict(test_label',double(fvt_test{i}'),model);
        end
        
        ac
        save(pred_results,'pred','ac','-v7.3');
        
        fprintf('Done!\n');
    case 'bow'
        %% Generate codebook
        cb_size=4000;
        redo_num=8;
        [Chog, Dhog, Ihog] = yael_kmeans(single(dtf_train{1}),cb_size,'redo', redo_num);
        [Chof, Dhof, Ihof] = yael_kmeans(single(dtf_train{2}),cb_size,'redo', redo_num);
        [Cmbhx, Dmbhx, Imbhx] = yael_kmeans(single(dtf_train{3}),cb_size,'redo', redo_num);
        [Cmbhy, Dmbhy, Imbhy] = yael_kmeans(single(dtf_train{4}),cb_size,'redo', redo_num);
        
        hog_cb=zeros(cb_size,size(dtf_train{1},2));
        hof_cb=zeros(cb_size,size(dtf_train{2},2));
        mbhx_cb=zeros(cb_size,size(dtf_train{3},2));
        mbhy_cb=zeros(cb_size,size(dtf_train{4},2));
        % Vector Quantization
        for i=1:size(dtf_train{1},2)
            % TBD %
        end
        
    otherwise
        error('Unsupported encoding method!');
end
