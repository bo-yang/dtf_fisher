function val = readIndex(obj, i, j)
%READINDEX Summary of this function goes here
%   Detailed explanation goes here

    linidx = obj.getIndex_(i, j);
    
    val = obj.matstore_(linidx);
end

