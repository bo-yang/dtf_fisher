function im = standardizeImage(im)
%STANDARDIZEIMAGE Summary of this function goes here
%   NOTE: all PASCAL VOC images are RGB and already size-normalized, so the
%   output of this function for them is always equivalent to:
%   single(rgb2gray(im))

%if ndims(im) == 3
%    im = rgb2gray(im);
%elseif ndims(im) ~= 2
%    error('Input image not valid');
%end

%im = single(im);

if ndims(im) == 3
    im = im2single(im);
elseif ndims(im) == 2
    im_new = cat(3,im,im);
    im_new = cat(3,im_new,im);
    im = im_new;
    im = im2single(im);
    clear im_new;
else
    error('Input image not valid');
end

if size(im,1) > 480, im = imresize(im, [480 NaN]) ; end

end

