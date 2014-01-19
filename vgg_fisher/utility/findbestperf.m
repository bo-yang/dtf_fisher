function [map, file] = findbestperf(pathloc, strmask)
%FINDBESTPERF Summary of this function goes here
%   Detailed explanation goes here

    if nargin == 0
        [files, path] = uigetfile('*.mat','Select result files',pwd,'MultiSelect','on');
        for i = 1:length(files)
            files{i} = fullfile(path, files{i});
        end
    if nargin == 1
        files = getFilenames(pathloc, {'mat'});
    elseif nargin == 2
        files = getFilenames(pathloc, {'mat'});
        for i = length(files):-1:1
            if isempty(strfind(files{i}, strmask))
                files(i) = [];
            end
        end
    else
        error('findbestperf requires between 0 and 2 arguments');
    end
    
    maxMap = 0;
    maxCFile = '';
    
    for i = 1:length(files)
        file = load(files{i});
        fmap = mean(file.results.res{1});
        if fmap > maxMap
            maxMap = fmap;
            maxCFile = files{i};
        end
    end
    
    map = maxMap;
    file = maxCFile;
end

