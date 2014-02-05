function [ fv ] = fisher_encode( feat, pca_coeff, gmm )
%FISHER_ENCODE compute fisher encoding by trained PCA coefficients and GMM
%parameters.

% L1 normalization & Sqare root
feat=sqrt(feat/norm(feat,1));

feat=pca_coeff*feat; % Apply PCA
fv=vl_fisher(double(feat),gmm.means,gmm.covar,gmm.prior,'Improved'); % Calc Fisher vector

end

