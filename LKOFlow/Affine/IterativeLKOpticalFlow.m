%d=IterativeLKOpticalFlow(img1, img2, roi, dInit)
%
% Computes the Lucas-Kanade optical flow using the iterative method.
% The function computes the flow, in each iteration the flow is refined.
% The iteration is stopped if the flow magnitude falls beneath 0.01 pixels
% or a maximum of 10 iterations occur.
%
% Inputs:
% img1 - The source image.
% img2 - The destination image.
% roi  - A region in img1 in which the flow should be computed.
% dInit - The initial estimate of the affine parameters
%
% Outpus:
% d - The computed affine transformation parameters
function d=IterativeLKOpticalFlow(img1, img2, roi, dInit)

K=10; % Number of iterations
STOP_THR = 0.01; % Stop if accuracy is better than 0.01 pixel

% Copy inflated region of interest our of image (we need border for
% derivative operation)
img1 = img1(roi(1)-1:roi(3)+1,roi(2)-1:roi(4)+1);

% Compute Ix and Iy, G and H
[Ht,G]=ComputeLKFlowParms(img1);

% Remove border
img1=img1(2:end-1, 2:end-1);

% Loop and perform iterations of optical flow computation
d=dInit;

k=1;

dc=inf;
while k<K && norm(dc)>STOP_THR

  % Use current translation vector to resample img2 and compute
  % temporal difference image
  It=img1-ResampleImg(img2, roi, d);

  % Find indexes which are not out of fov
  I = ~isnan(It(:));
  
  % Compute right side of optical flow equation: Gd=b
  b=Ht(:,I)*It(I);

  % Solve for current affine displacement parameters
  dc = G\b;
  
  % Add current displacement to d (This is actually concatinating the two
  % affine matrixes)
  d=d+reshape(dc, 2, 3)*(eye(3)+[d;0 0 0]);
  
  k=k+1;

end

dummy=1;