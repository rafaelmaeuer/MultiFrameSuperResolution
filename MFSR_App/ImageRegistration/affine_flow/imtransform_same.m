function im = imtransform_same(im, t, interp)
%IMTRANSFORM_SAME transforms an image into the original coordinates
%   IM = IMTRANSFORM_SAME(IM, T) transforms image IM using transform T,
%   produced e.g. with MAKETFORM, by calling IMTRANSFORM. Unlike the
%   default for IMTRANSFORM there is no change of coordinates - the origin
%   of the both the original and the new image is at row 0, column 0.
%
%   IM = IMTRANSFORM_SAME(IM, T, INTERP) passes INTERP as the third
%   argument to IMTRANSFORM.
%
%   See also: MAKETFORM, IMTRANSFORM

% Copyright David Young 2010

if nargin < 3
    im = imtransform(im, t, 'size', size(im), ...
        'xdata', [1 size(im,2)], 'ydata', [1 size(im, 1)]);
else
    im = imtransform(im, t, interp, 'size', size(im), ...
        'xdata', [1 size(im,2)], 'ydata', [1 size(im, 1)]);
end

end