function [ap, rec, prec] = computeAP(out, gt, poslabel, neglabel, draw)

% Given a list of confidences (in out) and the corresponding gt, compute AP
% Code is based off VOCevalcls

if ~exist('draw', 'var') || isempty(draw)
    draw = false;
end

if ~exist('poslabel', 'var') || isempty(poslabel)
    poslabel = 1;
end

if ~exist('neglabel', 'var') || isempty(neglabel)
    neglabel = -1;
end

% compute precision/recall

[~, si] = sort(-out);
tp = gt(si)==poslabel;
fp = gt(si)==neglabel;

fp  = cumsum(fp);
tp  = cumsum(tp);
rec = tp/sum(gt>0);
prec= tp./(fp+tp);

% compute average precision

ap = 0;
for t = 0:0.1:1
    p=max(prec(rec>=t));
    if isempty(p)
        p = 0;
    end
    ap = ap+p/11;
end

if draw
    % plot precision/recall
    plot(rec,prec,'-');
    grid;
    xlabel 'recall'
    ylabel 'precision'
    title(sprintf('AP = %.3f',ap));
end
