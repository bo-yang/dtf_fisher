function mat = convToFullMatrix(obj)
%CONVTOFULLMATRIX Summary of this function goes here
%   Detailed explanation goes here

    mat = zeros(obj.N);
    
    for i = 1:obj.N
        for j = 1:obj.N
            if i == j, continue; end
            if i < j, continue; end
            val = obj.readIndex(i, j);
            mat(i,j) = val;
            mat(j,i) = val;
        end
    end
end

