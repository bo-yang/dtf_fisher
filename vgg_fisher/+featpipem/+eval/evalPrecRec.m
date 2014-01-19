function AP = evalPrecRec(imdb, scoremat, testset, dataset)
%EVALPRECREC Evaluate a dataset using average precision
%   If the last parameter 'dataset' is not specified or is 'VOC2010' or
%   greater, the VOC2010+ precision-recall calculation methodology is used
%   (calculating precision continuously). If dataset is 'VOC2007-VOC2009',
%   the VOC2007-VOC2009 precision-recall calculation methodology is used
%   (calculating precision at fixed 0.1 intervals)

if nargin < 4 || isempty(dataset)
    dataset = '';
end

AP = zeros(1,size(scoremat,1));
%AP = struct;

fprintf('Retrieving ground truth from IMDB...\n');
gt = featpipem.utility.getImdbGT(imdb, {testset}, 'outputSignedLabels', true);

for ci = 1:size(scoremat,1)    
    % get ground truth
    gt_cls = gt{ci};
    if length(gt_cls) ~= size(scoremat,2)
        error('Mismatch between size of scoremat and number of images in testset(s)');
    end
    
    % get indices of current class sorted in descending order of confidence
    fprintf('Computing precision/recall...\n');
    [sortidx sortidx] = sort(scoremat(ci,:),'descend'); %#ok<ASGLU>
    tp = gt_cls(sortidx)>0;
    fp = gt_cls(sortidx)<0;
    
    fp = cumsum(fp);
    tp = cumsum(tp);
    rec = tp/sum(gt_cls>0);
    prec = tp./(fp+tp);
    
    fprintf('Computing APs...\n');
    AP(ci) = featpipem.eval.VOCdevkit.VOCap(dataset, rec, prec);    
    %AP(ci).AP = featpipem.eval.VOCdevkit.VOCap(dataset, rec, prec);
    %AP(ci).prec10 = prec(10);
    %AP(ci).prec50 = prec(50);
    
    fprintf('DONE\n');
    
end

end

