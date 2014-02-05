function [ labels, file_list ] = load_dtf_data( params, file_type )
%LOAD_DTF_DATA Load DTF data
%   Assume that the dtf files have similar hierachy with UCF 101
%   Inputs:
%       dtf_dir - directory of gzipped DTF descriptors
%       tt_list_dir - directory of train/test file list
%       file_pattern - pattern of train/test list, such as 'train*' or 'test*'
%       num_feats - number of features extracted from each dtf file
%
%   Outputs:
%       L -labels for DTF features

switch file_type
    case 'train'
        tt_list_dir=params.train_list_dir;
        reg_pattern='train*';
    case 'test'
        tt_list_dir=params.test_list_dir;
        reg_pattern='test*';
    otherwise
        error('Unknown file pattern!');
end

% extract training/test list
tt_list=[]; % train/test files
labels=[]; % labels
tlists=dir(fullfile(tt_list_dir,reg_pattern));
for i=1:length(tlists)
    fid=fopen(fullfile(tt_list_dir,tlists(i).name));
    tmp=textscan(fid,'%s%d');
    tt_list=[tt_list;tmp{1}];
    labels=[labels;tmp{2}];
end

file_list=cell(length(tt_list),1);
for i=1:length(tt_list)
    clip_name=regexprep(tt_list{i},'\.avi$',''); % get video clip name
    clip_name=regexprep(clip_name,'.*/','');
    file_list{i}=[clip_name,'.txt'];
end

% extract DTF data and set labels
%try
%	matlabpool close;
%catch exception
%end
%matlabpool open 4
parfor i=1:length(tt_list)
	action=regexprep(tt_list{i},'/v_(\w*)\.avi','');
    act_dir=fullfile(params.dtf_dir,action);
    clip_name=regexprep(tt_list{i},'\.avi$',''); % get video clip name
    clip_name=regexprep(clip_name,'.*/','');
    dtf_file=fullfile(act_dir,[clip_name,'.dtf.gz']);
	saved_feat_file=sprintf('%s.txt',clip_name);
	
	switch file_type
    case 'train'
        HOG_file=fullfile(params.HOG_train_data,saved_feat_file);
        HOF_file=fullfile(params.HOF_train_data,saved_feat_file);
        MBHx_file=fullfile(params.MBHx_train_data,saved_feat_file);
        MBHy_file=fullfile(params.MBHy_train_data,saved_feat_file);
    case 'test'
        HOG_file=fullfile(params.HOG_test_data,saved_feat_file);
        HOF_file=fullfile(params.HOF_test_data,saved_feat_file);
        MBHx_file=fullfile(params.MBHx_test_data,saved_feat_file);
        MBHy_file=fullfile(params.MBHy_test_data,saved_feat_file);
    otherwise
        error('Unknown file pattern!');
	end
    
    [HOG,HOF,MBHx,MBHy]=extract_dtf_feats(dtf_file,params);
    % To save memory, write DTF features into txt files instead of stored in memory
    dlmwrite(HOG_file,HOG,'delimiter',' ');
    dlmwrite(HOF_file,HOF,'delimiter',' ');
    dlmwrite(MBHx_file,MBHx,'delimiter',' ');
    dlmwrite(MBHy_file,MBHy,'delimiter',' ');
end
%matlabpool close

end

