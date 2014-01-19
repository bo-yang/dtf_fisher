function linidx = getIndex_(obj, i, j)
%GETINDEX Summary of this function goes here
%   Detailed explanation goes here

    if i == j
        error('No diagonal specified in lower triangular matrix');
    end
    % ensure i is always greater than j
    if i < j
        tmp = i;
        i = j;
        j = tmp;
    end
    
    % real_i and real_j refer to index in LT matrix
    % i.e.
    % - - -        * * *
    % * - -    ==> * *
    % * * -        *
    % * * *
    real_i = i-j;
    real_j = j;
    
    % calculate index in collapsed vector
    linidx = 0;
    coladd = obj.N-1;
    for col = 1:(real_j-1)
        linidx = linidx + coladd;
        coladd = coladd - 1;
    end
    linidx = linidx+real_i;
end

