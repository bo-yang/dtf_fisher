function [ pca_coeff, gmm, file_list, labels ] = pretrain(params, t_type)
%PRETRAIN  Subsample DTF features, calculate PCA coefficients and train GMM model
%   Inputs:
%       params - structure of parameters
%       t_type - 'train' or 'test'
%
%   Outputs:
%       pca_coeff - PCA coefficients for each DTF feature
%       gmm - GMM params for each DTF feature

% To construct Fisher vector, fist to estimate the parameters of GMM.
% To estimate GMM params, subsampling vectors from DTF descriptors.

gmm_params.cluster_count=params.K;
gmm_params.maxcomps=gmm_params.cluster_count/4;
gmm_params.GMM_init= 'kmeans';
gmm_params.pnorm = single(2);    % L2 normalization, 0 to disable
gmm_params.subbin_norm_type = 'l2';
gmm_params.norm_type = 'l2';
gmm_params.post_norm_type = 'none';
gmm_params.pool_type = 'sum';
gmm_params.quad_divs = 2;
gmm_params.horiz_divs = 3;
gmm_params.kermap = 'hellinger';

dtf_feat_num=length(params.feat_list);
pca_coeff=cell(dtf_feat_num,1); % PCA coeficients
gmm=cell(dtf_feat_num,1);   % GMM parameters

fprintf('\nPreprocessing features for %s ...\n', t_type);

switch t_type
	case 'train'
		feat_sample_file=params.train_sample_data;
	case 'test'
		feat_sample_file=params.test_sample_data;
	otherwise
		error('Only TRAIN or TEST are valid types!');
end

fprintf('Subsampling DTF features ...\n');
if ~exist(feat_sample_file,'file')
	[feats, file_list, labels]=subsample(params,t_type);
	save(feat_sample_file, 'feats','labels','file_list','-v7.3');
else
	load(feat_sample_file);
end

for i=1:dtf_feat_num
	feat=feats{i};
	
	% L1 normalization & Sqare root
	feat=sqrt(feat/norm(feat,1));
	
	% Do PCA on train/test data to half-size original descriptors
	fprintf('Doing PCA ...\n');
	pca_coeff{i} = princomp(feat');
	pca_coeff{i} = pca_coeff{i}(:, 1:floor(size(feat,1)/2))';
	% dimensionality reduction
	feat = pca_coeff{i} * feat;

	fprintf('Training Guassian Mixture Model ...\n');
	%[gmm{i}.means, gmm{i}.covar, gmm{i}.prior] = vl_gmm(feat, params.K, 'MaxNumIterations', 300);
	gmm{i}=gmm_gen_codebook(feat,gmm_params);
end

clear feats

end

