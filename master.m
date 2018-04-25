% Script to test DTF features and Fisher Vector

run('vlfeat-0.9.17/toolbox/vl_setup')

path=pwd;
addpath(path);
%addpath(fullfile(path,'yael'));
%addpath(fullfile(path,'yael/matlab'));
addpath(fullfile(path,'libsvm'));
addpath(fullfile(path,'libsvm/matlab'));
addpath(fullfile(path,'vgg_fisher'));
addpath(fullfile(path,'vgg_fisher/lib/gmm-fisher/matlab'));

%% Set Parameters
params.K=256;   % num of GMMs
params.DTF_subsample_num=1000; % Subsampling number of DTF features per video clip

params.encoding='fisher'; % encoding type: 'fisher' - fisher vector; 'bow' - bag-of-words
params.feat_list={'HOG','HOF','MBHx','MBHy'}; % all features involved in this test
feat_len={96,108,96,96}; % length of features
params.feat_len_map=containers.Map(params.feat_list, feat_len);
params.feat_start=41; % start position of DTF features

params.dtf_dir=fullfile(path,'UCF101_DTF');
params.train_list_dir=fullfile(path, 'ucfTrainTestlist'); 
params.test_list_dir=fullfile(path, 'ucfTrainTestlist'); 

params.train_data_info=fullfile(path,'data','UCF101_traindata_info.mat');
params.test_data_info=fullfile(path,'data','UCF101_testdata_info.mat');

% Files to store subsampled features
params.train_sample_data=fullfile(path,'data',sprintf('UCF101_train_data_sample%d_gmm%d.mat',params.DTF_subsample_num,params.K));
params.test_sample_data=fullfile(path,'data',sprintf('UCF101_test_data_sample%d_gmm%d.mat',params.DTF_subsample_num,params.K));

% Files to store Fisher vectors
params.fv_train_file=fullfile(path,'data',sprintf('pca_gmm_data_train_sample%d_gmm%d.mat',params.DTF_subsample_num,params.K));
params.fv_test_file=fullfile(path,'data',sprintf('pca_gmm_data_test_sample%d_gmm%d.mat',params.DTF_subsample_num,params.K));

pred_results=fullfile(path,'data',sprintf('pred_results_sample%d_gmm%d.mat',params.DTF_subsample_num,params.K));

try
	matlabpool close;
catch exception
end
matlabpool open 4

switch params.encoding
      case 'fisher'        
        %% Construct Fisher Vectors for training
        % Subsample DTF features, calculate PCA coefficients and train GMM model
        if ~exist(params.fv_train_file,'file') && ~exist(params.fv_test_file,'file')
            % Train and test share the same codebook?
            [ pca_coeff, gmm, all_train_files, all_train_labels, all_test_files, all_test_labels ] = pretrain(params);
            save(params.fv_train_file,'pca_coeff','gmm','all_train_files','all_train_labels','all_test_files', 'all_test_labels', '-v7.3');
        else
            load(params.fv_train_file); 
        end
        
        % Load trainning videos, compute Fisher vectors and train SVM model
        svm_option='-t 0 -s 0 -q -w0 0.5 -w1 0.5 -c 100 -b 1'; % temporarily not to use linear SVM
        uniq_labels=unique(all_train_labels);
        pred=zeros(length(all_test_files),1);
        acc=zeros(numel(uniq_labels),1);
        %pred=[];
        
        for i=1:numel(uniq_labels)
            % Process training files
            pos_idx=(all_train_labels == i);
            pos_files=all_train_files(pos_idx); % positive training files
            fvt_pos_train=compute_fisher(params, pca_coeff, gmm, pos_files);
            
            neg_idx=(all_train_labels ~= i);
            neg_files=all_train_files(neg_idx);
            sample_idx=randperm(length(neg_files),length(pos_files));
            sample_neg_files=neg_files(sample_idx);
            fvt_neg_train=compute_fisher(params, pca_coeff, gmm, sample_neg_files);
            
            % Train SVM model
            tmp_train_labels=[ones(1,size(fvt_pos_train,2)) -1*ones(1,size(fvt_neg_train,2))];
            model=svmtrain(double(tmp_train_labels)', double([fvt_pos_train fvt_neg_train])', svm_option);
            
            % Process test files
            pos_idx=(all_test_labels == i);
            pos_files=all_test_files(pos_idx); % positive training files
            fvt_pos_test=compute_fisher(params, pca_coeff, gmm, pos_files);
            
            neg_idx=(all_test_labels ~= i);
            neg_files=all_test_files(neg_idx);
            sample_idx=randperm(length(neg_files),length(pos_files));
            sample_neg_files=neg_files(sample_idx);
            fvt_neg_test=compute_fisher(params, pca_coeff, gmm, sample_neg_files);
            
            % SVM prediction
            tmp_test_labels=[ones(1,size(fvt_pos_test,2)) -1*ones(1,size(fvt_neg_test,2))];
            [pred_labels,accuracy,prob_estimates] = svmpredict(double(tmp_test_labels)', double([fvt_pos_test fvt_neg_test])', model, '-b 1');
            acc(i) = sum(pred_labels(1:size(fvt_pos_test,2)) == 1) ./ size(fvt_pos_test,2);    %# accuracy
            pred(pos_idx) = (pred_labels(1:size(fvt_pos_test,2))==1)*i;
            %pred = [pred; (pred_labels(1:size(fvt_pos_test,2))==1)*i];
            
			save_file=sprintf('fisher_vectors_action%d_sample%d_gmm%d.mat',i,params.DTF_subsample_num,params.K);
			save(fullfile(path,'data',save_file),'fvt_pos_train','fvt_neg_train', 'fvt_pos_test', 'fvt_neg_test', '-v7.3');
        end

        %fprintf('Mean accuracy: %f.\n', mean(acc)); % display mean accuracy 
		fprintf('\nMean accuracy: %f.\n', sum(pred==all_test_labels)/numel(all_test_labels)); % display mean accuracy 
        save(pred_results,'pred','acc','-v7.3');
        
        
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

try
	matlabpool close;
catch exception
end

fprintf('Done!\n');
