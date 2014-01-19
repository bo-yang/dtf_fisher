function writeIndex(obj, i, j, val)
%WRITEINDEX Summary of this function goes here
%   Detailed explanation goes here

    linidx = obj.getIndex_(i, j);
    
    obj.matstore_(linidx) = val;
end

