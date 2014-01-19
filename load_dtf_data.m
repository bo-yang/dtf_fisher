function [ X, L ] = load_dtf_data( dtf_dir, tt_list_dir, file_pattern, encoding, num_feats )
%LOAD_DTF_DATA Load DTF data
%   Assume that the dtf files have similar hierachy with UCF 101
%   Inputs:
%       dtf_dir - directory of gzipped DTF descriptors
%       tt_list_dir - directory of train/test file list
%       file_pattern - pattern of train/test list, such as 'train*' or 'test*'
%       num_feats - number of features extracted from each dtf file
%
%   Outputs:
%       X - DTF features
%       L -labels for DTF features

if ~exist('num_feats','var')
    num_feats=100; % default number of features extracted from each video
end

% extract training/test list
tt_list=[]; % train/test files
ll_list=[]; % labels
tlists=dir(fullfile(tt_list_dir,file_pattern));
for i=1:length(tlists)
    fid=fopen(fullfile(tt_list_dir,tlists(i).name));
    tmp=textscan(fid,'%s%d');
    tt_list=[tt_list;tmp{1}];
    ll_list=[ll_list;tmp{2}];
end

X=cell(4,1);
%L=zeros(1,length(tt_list)*num_feats);
L=[];

% extract DTF data and set labels
for i=1:length(tt_list)
    action=regexprep(tt_list{i},'/v_(\w*)\.avi','');
    act_dir=fullfile(dtf_dir,action);
    clip_name=regexprep(tt_list{i},'\.avi$',''); % get video clip name
    clip_name=regexprep(clip_name,'.*/','');
    dtf_file=fullfile(act_dir,[clip_name,'.dtf.gz']);
    
    if ~exist(dtf_file,'file')
        warning('File %s does not exist! Skip now...',dtf_file);
        continue;
    else
        tmpfile=dir(dtf_file);
        if tmpfile.bytes < 1024
            warning('File %s is too small! Skip now...',dtf_file);
            continue;
        end
    end
    
    [HOG,HOF,MBHx,MBHy]=extract_dtf_feats(dtf_file,num_feats);
    X{1}=[X{1} HOG];
    X{2}=[X{2} HOF];
    X{3}=[X{3} MBHx];
    X{4}=[X{4} MBHy];
    
    switch encoding
        case 'bow'
            %L((i-1)*num_feats+1:i*num_feats)=ll_list(i);
            label=zeros(1,num_feats)+double(ll_list(i));
            L=[L label];
        case 'fisher'
            L=[L ll_list(i)];
        otherwise
            error('unsupported encoding method!');
    end
end

end

