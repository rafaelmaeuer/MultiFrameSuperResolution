% Implements the robust super-resolution method. This funtion uses the
% steepest descent method to minimize the SR cost function which includes
% two terms. The "energy" term, which is the L1 norm of the residual error
% between the HR image and the LR image sequence. The "regulerization" term
% which induces piecewise smoothness on the HR image using the bilteral
% filter.
%
% Inputs:
% LR - A sequence of low resolution images
% D  - The tranlational motion for each LR frame
% HR - An initial guess for the HR image.
% resFactor - The resolution increment factor
% Hpsf - The PSF function (common to all frames and space invariant)
% props - property structure used to control the algorithm parameters
%
% Outpus:
% The estimated HR image
function HR=RobustSR(LR, D, HR, resFactor, Hpsf, props)

% Loop and improve HR in steepest descent direction
iter = 1;

h=waitbar(0, 'Estimating high-resolution image');

while iter<props.maxIter
  
  waitbar(iter/props.maxIter);
  
  % Compute gradient of the energy part of the cost function
  Gback = GradientBackProject(HR, LR, D, Hpsf, resFactor);
  
  % Compute the gradient of the bilateral filter part of the cost function
  Greg = GradientRegulization(HR, props.P, props.alpha);

  % Perform a single SD step
  HR = HR - props.beta.*(Gback + props.lambda.* Greg);
  
  iter=iter+1;

end

close(h);