% Perform Fisher Vector encoding
function code=fisher_encode_vgg(feats,pca_coeff, codebook, params, weights)

% L1 normalization & Sqare root
feats=sqrt(feats/norm(feats,1));
feats=pca_coeff*feats; % Apply PCA

 % Initialize encoder
 if (nargin < 4) || isempty(params)
    cpp_handle = mexFisherEncodeHelperSP('init', codebook);
 else
    cpp_handle = mexFisherEncodeHelperSP('init', codebook, params);
 end
 
 % Encode
 if (nargin < 5) || isempty(weights)
     code = mexFisherEncodeHelperSP('encode', cpp_handle, single(feats));
 else
     code = mexFisherEncodeHelperSP('encode', cpp_handle, single(feats), weights);
 end

  % destructor
  mexFisherEncodeHelperSP('clear', cpp_handle);


%         % get FK dimensionality
%         function dim = getdim(this)
%             dim = mexFisherEncodeHelperSP('getdim', this.cpp_handle);
%         end
%     end
%     
% end


end
