function [ HOG,HOF,MBHx,MBHy ] = extract_dtf_feats( dtf_file, num_feats)
%EXTRACT_DTF_FEATS extract DTF features.
%   The first 10 elements for each line in dtf_file are information about the trajectory.
%   The trajectory info(default 30 dimensions) should also be discarded.
%
%   Subsampling: 
%       randomly choose 100 descriptors from each video clip(dtf file)
%

if ~exist('num_feats','var')
    num_feats=100;
end

unzip_cmd=sprintf('gunzip %s',dtf_file); % Suppose dtf files have suffix .gz
system(unzip_cmd); 
unzip_dtf_file=regexprep(dtf_file,'\.gz$',''); % remove suffix .gz
x=load(unzip_dtf_file);
zip_cmd=sprintf('gzip -f %s',unzip_dtf_file); % Suppose dtf files have suffix .gz
system(zip_cmd); 

feat_start=41;
hog_range=feat_start:feat_start+96-1;
hof_range=hog_range(end)+1:hog_range(end)+108;
mbhx_range=hof_range(end)+1:hof_range(end)+96;
mbhy_range=mbhx_range(end)+1:mbhx_range(end)+96;

HOG=zeros(numel(hog_range),num_feats);
HOF=zeros(numel(hof_range),num_feats);
MBHx=zeros(numel(mbhx_range),num_feats);
MBHy=zeros(numel(mbhy_range),num_feats);

if size(x,1)<num_feats
    warning('File %s has fewer rows than specified num_feats!', dtf_file);
    % Replicate x so as to make x large enough for subsampling
    x=repmat(x,ceil(num_feats/size(x,1)),1);
end

idx=randperm(size(x,1),num_feats); % randomly subsampling
%idx=floor(linspace(1,size(x,1),num_feats)); % linearly subsampling
HOG=x(idx,hog_range)';
HOF=x(idx,hof_range)';
MBHx=x(idx,mbhx_range)';
MBHy=x(idx,mbhy_range)';

end

