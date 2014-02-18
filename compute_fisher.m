function [ fvt ] = compute_fisher( params, pca_coeff, gmm, file_list )
%COMPUTE_FISHER Computer postitive/negative Fisher Vectors for each class

fvt=[];

fisher_params.alpha = single(0.5);  % power notmalization, 1 to disable
fisher_params.grad_weights = false; % soft BOW
fisher_params.grad_means = true;    % 1st order
fisher_params.grad_variances = true;    % 2nd order

% Make sure all files are compressed
check_gzip_cmd=sprintf('sh ./gzip_dtf_files -i %s > /dev/null 2>&1',params.dtf_dir);
system(check_gzip_cmd);

%feat_idx=find(strcmp(params.feat_list,feat_type)); % find the index of feature in feat_list

% TODO:
% Should the following code be executed in parallel?
% The label of each FV is not an issue, because all of them could be
% assigned 1(positive files) or -1(negative files).
% What about the concatenation of Fisher vectors?
for j=1:length(file_list)
	action=regexprep(file_list{j},'/v_(\w*)\.avi','');
	act_dir=fullfile(params.dtf_dir,action);
	clip_name=regexprep(file_list{j},'\.avi$',''); % get video clip name
	clip_name=regexprep(clip_name,'.*/','');
	file=fullfile(act_dir,[clip_name,'.dtf.gz']);
	
	[HOG,HOF,MBHx,MBHy]=extract_dtf_feats(file, params, -1);
	
	%fv_hog=fisher_encode(HOG,pca_coeff{1},gmm{1});
	%fv_hof=fisher_encode(HOF,pca_coeff{2},gmm{2});
	%fv_MBHx=fisher_encode(MBHx,pca_coeff{3},gmm{3});
	%fv_MBHy=fisher_encode(MBHy,pca_coeff{4},gmm{4});
	fv_hog=fisher_encode_vgg(HOG,pca_coeff{1},gmm{1},fisher_params);
	fv_hof=fisher_encode_vgg(HOF,pca_coeff{2},gmm{2},fisher_params);
	fv_MBHx=fisher_encode_vgg(MBHx,pca_coeff{3},gmm{3},fisher_params);
	fv_MBHy=fisher_encode_vgg(MBHy,pca_coeff{4},gmm{4},fisher_params);
	
	fv=[fv_hog;fv_hof;fv_MBHx;fv_MBHy];
	fvt=[fvt fv]; % concatenate all features together
end

% power normalization
fvt = sign(fvt) .* sqrt(abs(fvt));
% L2 normalization
fvt = double(yael_fvecs_normalize(single(fvt)));

fvt(find(isnan(fvt))) = 123456;

end

