% Implements the fast and robust super-resolution method. This funtion
% first compute an estimation of the blurred HR image, using the median and
% shift method. It then uses the bilateral filter as a regulating term
% for the debluring and interpolation step.
%
% Inputs:
% LR - A sequence of low resolution images
% D  - The tranlational motion for each LR frame
% resFactor - The resolution increment factor
% Hpsf - The PSF function (common to all frames and space invariant)
% props - property structure used to control the algorithm parameters
%
% Outputs:
% The estimated HR image
function HR=FastRobustSR(LR, D, resFactor, Hpsf, props)

% Compute initial estimate of blurred HR by the means of MedianAndShift
[Z, A]=MedianAndShift(LR, D, [(size(LR,1)+1)*resFactor-1 (size(LR,2)+1)*resFactor-1], resFactor);

% Deblur the HR image and regulate using bilatural filter

% Loop and improve HR in steepest descent direction
HR = Z;
iter = 1;

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