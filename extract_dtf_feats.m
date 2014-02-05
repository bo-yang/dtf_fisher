function [ HOG,HOF,MBHx,MBHy ] = extract_dtf_feats( dtf_file, params, num_feats )
%EXTRACT_DTF_FEATS extract DTF features.
%   The first 10 elements for each line in dtf_file are information about the trajectory.
%   The trajectory info(default 30 dimensions) should also be discarded.
%	
%   Subsampling:
%       randomly choose 100 descriptors from each video clip(dtf file)
%		To use all the DTF fatures, set num_feats to a negative number.
%

HOG=zeros(params.feat_len_map('HOG'),1);
HOF=zeros(params.feat_len_map('HOF'),1);
MBHx=zeros(params.feat_len_map('MBHx'),1);
MBHy=zeros(params.feat_len_map('MBHy'),1);

if ~exist(dtf_file,'file')
	warning('File %s does not exist! Skip now...',dtf_file);
	return;
else
	tmpfile=dir(dtf_file);
	if tmpfile.bytes < 1024
		warning('File %s is too small! Skip now...',dtf_file);
		return;
	end
end

unzip_cmd=sprintf('gunzip %s',dtf_file); % Suppose dtf files have suffix .gz
system(unzip_cmd);
unzip_dtf_file=regexprep(dtf_file,'\.gz$',''); % remove suffix .gz
x=load(unzip_dtf_file);
zip_cmd=sprintf('gzip -f %s',unzip_dtf_file); % Suppose dtf files have suffix .gz
system(zip_cmd);

hog_range=params.feat_start:params.feat_start+params.feat_len_map('HOG')-1;
hof_range=hog_range(end)+1:hog_range(end)+params.feat_len_map('HOF');
mbhx_range=hof_range(end)+1:hof_range(end)+params.feat_len_map('MBHx');
mbhy_range=mbhx_range(end)+1:mbhx_range(end)+params.feat_len_map('MBHy');

if num_feats<0 % To use all the DTF fatures, set num_feats to a negative number
	num_feats=size(x,1);
end

if size(x,1)<=num_feats
	idx=1:size(x,1); % randomly subsampling
else
	idx=randperm(size(x,1),num_feats); % randomly subsampling
	%idx=floor(linspace(1,size(x,1),num_feats)); % linearly subsampling
end
HOG=x(idx,hog_range)';
HOF=x(idx,hof_range)';
MBHx=x(idx,mbhx_range)';
MBHy=x(idx,mbhy_range)';

end

