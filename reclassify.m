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

%% SVM classification
load(params.fv_train_file);
%load(params.fv_test_file);
uniq_labels=unique(all_train_labels);
pred=zeros(length(all_test_files),1);
acc=zeros(numel(uniq_labels),1);
svm_option='-t 0 -s 0 -q -w0 0.5 -w1 0.5 -c 0.1 -b 1'; % temporarily not to use linear SVM

for i=1:numel(uniq_labels)
	save_file=sprintf('fisher_vectors_action%d_sample%d_gmm%d.mat',i,params.DTF_subsample_num,params.K);
	load(fullfile(path,'data',save_file));
	
	pos_idx=(all_test_labels == i);
	tmp_train_labels=[ones(1,size(fvt_pos_train,2)) -1*ones(1,size(fvt_neg_train,2))];
	model=svmtrain(double(tmp_train_labels)', double([fvt_pos_train fvt_neg_train])', svm_option);
	
	tmp_test_labels=[ones(1,size(fvt_pos_test,2)) -1*ones(1,size(fvt_neg_test,2))];
	[pred_labels,accuracy,prob_estimates] = svmpredict(double(tmp_test_labels)', double([fvt_pos_test fvt_neg_test])', model, '-b 1');
	acc(i) = sum(pred_labels(1:size(fvt_pos_test,2)) == 1) ./ size(fvt_pos_test,2);    %# accuracy
	pred(pos_idx) = (pred_labels(1:size(fvt_pos_test,2))==1)*i;
end

fprintf('\nMean accuracy: %f.\n', sum(pred==all_test_labels)/numel(all_test_labels)); % display mean accuracy
