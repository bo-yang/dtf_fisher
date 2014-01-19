% Removes DC component from image patches
% Data given as a matrix where each patch is one column vectors
% That is, the patches are vectorized.

function Y=removeDC(X);

% Subtract local mean gray-scale value from each patch in X to give output Y

Y = X-ones(size(X,1),1)*mean(X);

return;
