% Perform Fisher Vector encoding
function code=fisher_encode(feats,codebook, params, weights)

 % Initialize encoder
 if (nargin < 3) || isempty(params)
    cpp_handle = mexFisherEncodeHelperSP('init', codebook);
 else
    cpp_handle = mexFisherEncodeHelperSP('init', codebook, params);
 end
 
 % Encode
 if (nargin < 4) || isempty(weights)
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
