function imdb = generateImdbFromFolders(pathloc)
%GENERATEIMDBFROMFOLDERS Summary of this function goes here
%   Detailed explanation goes here

    imdb.images.id = [];
    imdb.images.set = uint8([]);
    imdb.images.name = {};
    imdb.images.class = [];
    imdb.images.size = zeros(2,0) ;
    imdb.dir = getTopLevelDir(pathloc);
    imdb.classes.name = [];
    imdb.classes.imageIds = {};
    imdb.sets = struct();

    sets = getDirectoriesAtPath(pathloc);
    for i = 1:length(sets)
        if strcmpi(sets{i},'train')
            if i ~= 1
                tmp = sets{1};
                sets{1} = sets{i};
                sets{i} = tmp;
            end
        end
    end
    for i = 1:length(sets)
        if strcmpi(sets{i},'val')
            if i ~= 2
                tmp = sets{2};
                sets{2} = sets{i};
                sets{i} = tmp;
            end
        end
    end
    
    
    % store set names in imdb
    for si = 1:length(sets)
        imdb.sets.(sets{si}) = uint8(si);
    end
    
    imid = 1;
    
    for si = 1:length(sets)
        % get list of classes from subdirs in set folder
        classes = getDirectoriesAtPath(fullfile(pathloc, sets{si}));
        % check classes are same for all sets
        if isempty(imdb.classes.name)
            imdb.classes.name = classes;
            imdb.classes.imageIds = cell(1, length(imdb.classes.name));
        else
            if ~isempty(setxor(classes, imdb.classes.name))
                error('Classes must be same for all datasets');
            end
        end
        % enter each directory in turn, adding images to the imdb
        for ci = 1:length(classes)
            ims = getImagesAtPath(fullfile(pathloc, sets{si}, classes{ci}));
            for ii = 1:length(ims)
                % load in image to get size
                im = imread(fullfile(pathloc, sets{si}, classes{ci}, ims{ii}));
                imsize = size(im)';
                imsize = imsize(1:2);
                % store image in imdb
                imdb.images.name = [imdb.images.name {fullfile(sets{si}, classes{ci}, ims{ii})}];
                imdb.images.id = [imdb.images.id imid];
                imdb.images.set = [imdb.images.set uint8(si)];
                imdb.images.class = [imdb.images.class ci];
                imdb.images.size = [imdb.images.size imsize];
                
                imdb.classes.imageIds{ci} = [imdb.classes.imageIds{ci} imid];
                
                imid = imid+1;
            end
        end
    end
end

function dirs = getDirectoriesAtPath(pathloc)
    dirlisting_all = dir(pathloc);
    % remove '.' and '..' items from directory listing
    dirlisting_all(arrayfun(@(x)(strcmp(x.name,'.')),dirlisting_all)) = [];
    dirlisting_all(arrayfun(@(x)(strcmp(x.name,'..')),dirlisting_all)) = []; 
    %remove all non-directory entries
    dirlisting = dirlisting_all(arrayfun(@(x)(x.isdir),dirlisting_all));
    dirs = arrayfun(@(x)(x.name),dirlisting,'UniformOutput',false);
end

function ims = getImagesAtPath(pathloc)
    dirlisting_all = dir(pathloc);
    % remove '.' and '..' items from directory listing
    dirlisting_all(arrayfun(@(x)(strcmp(x.name,'.')),dirlisting_all)) = [];
    dirlisting_all(arrayfun(@(x)(strcmp(x.name,'..')),dirlisting_all)) = []; 
    %remove all non-image entries
    extstr = 'jpg|jpeg|gif|png|bmp';
    dirlisting = arrayfun(@(x)~isempty(regexpi(x.name,extstr)),dirlisting_all);
    ims = arrayfun(@(x)(x.name),dirlisting_all(dirlisting),'UniformOutput',false);
end

function dirname = getTopLevelDir(pathloc)
    filesepids = strfind(pathloc,filesep);
    if filesepids(end) == length(pathloc)
        pathloc = pathloc(1:end-1);
        filesepids = filesepids(1:end-1);
    end
    dirname = pathloc(filesepids(end)+1:end);
end

