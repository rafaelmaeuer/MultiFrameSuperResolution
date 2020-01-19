% Implements a simple cubic-spline interpolation of a single image. This
% image is then deblured using the same method as in the Fast and Robust
% method.
%
% Inputs:
% LR - A sequence of low resolution images
% resFactor - The resolution increment factor
% Hpsf - The PSF function (common to all frames and space invariant)
% props - property structure used to control the algorithm parameters
%
% Outpus:
% The estimated HR image
function HR=SplineSRInterp(LR, resFactor, Hpsf, props)

% Initialize guess as interpolated version of LR
[X,Y]=meshgrid(0:resFactor:(size(LR,2)-1)*resFactor, 0:resFactor:(size(LR,1)-1)*resFactor);
[XI,YI]=meshgrid(resFactor+1:(size(LR,2)-2)*resFactor-1, resFactor+1:(size(LR,1)-2)*resFactor-1);

Z=interp2(X, Y, squeeze(LR(:,:,1)), XI, YI, '*spline');

% Deblur the HR image and regulate using bilatural filter

% Loop and improve HR in steepest descent direction
HR = Z;
iter = 1;
A = ones(size(HR));

h=waitbar(0, 'Estimating high-resolution image');

while iter<props.maxIter
  
  waitbar(iter/props.maxIter);
  
  % Compute gradient of the energy part of the cost function
  Gback = FastGradientBackProject(HR, Z, A, Hpsf);

  % Compute the gradient of the bilateral filter part of the cost function
  Greg = GradientRegulization(HR, props.P, props.alpha);

  % Perform a single SD step
  HR = HR - props.beta.*(Gback + props.lambda.* Greg);
  
  iter=iter+1;

end

close(h);