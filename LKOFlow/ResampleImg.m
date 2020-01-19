% img=ResampleImg(img, roi, d) 
%
% Resamples the image at the given ROI and displacement
function img=ResampleImg(img, roi, d)

% Compute coordinates for resampling
[X,Y]=meshgrid((roi(2):roi(4))-d(1), (roi(1):roi(3))-d(2));

% Resample the image
img=interp2(img, X, Y);
