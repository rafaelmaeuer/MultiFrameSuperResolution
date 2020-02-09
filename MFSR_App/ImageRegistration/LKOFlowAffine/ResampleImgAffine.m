% img=ResampleImgAffine(img, roi, d) 
%
% Resamples the image at the given ROI and displacement
function img=ResampleImgAffine(img, roi, d)

% Compute coordinates for resampling
[X,Y]=meshgrid(roi(2):roi(4), roi(1):roi(3));

% Add affine displacement to coordinates
temp = [X(:)';Y(:)' ; ones(1,length(X(:)))];
D=d*temp;

% Update coorindates with displacement
X = X-reshape(D(1,:), size(X));
Y = Y-reshape(D(2,:), size(X));

% Resample the image
img=interp2(img, X, Y);
