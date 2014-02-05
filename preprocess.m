function [feats,train_labels,feat_nums] = preprocess(file_label_map,params,feat_type,t_type)

feat_len=params.feat_len_map(feat_type);
switch t_type
    case 'train'
        feat_dir=params.train_file_map(feat_type);
    case 'test'
        feat_dir=params.test_file_map(feat_type);
    otherwise
        error('Only TRAIN or TEST are valid types!');
end

feats=[];
train_labels=[];
feat_nums=[];
% Load features
all_files=keys(file_label_map);
try
	matlabpool close;
catch exception
end
matlabpool open 4
parfor i=1:length(all_files)
    file=fullfile(feat_dir,all_files{i});
    tmpfile=dir(file);
    try
        if tmpfile.bytes == 0
            file
            error('File contains no data!');
        end
    catch exception
        file
    end
    
    feat=dlmread(file);
    if size(feat,1) ~= feat_len
        file %%% TEST ONLY %%%
        size(feat) %%% TEST ONLY %%%
        feat=reshape(feat, floor(size(feat,1)*size(feat,2)/feat_len), feat_len)';
    end
    num_feats=ceil(size(feat,2)*params.DTF_subsample_rate);
    %idx=randperm(size(x,1),num_feats); % randomly subsampling
    idx=floor(linspace(1,size(feat,2),num_feats)); % linearly subsampling
    feats=[feats feat(:,idx)];
    train_labels=[train_labels file_label_map(all_files{i})];
    feat_nums=[feat_nums num_feats];
end
matlabpool close

end