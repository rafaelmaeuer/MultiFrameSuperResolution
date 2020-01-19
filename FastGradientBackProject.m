% Computes the gradient backprojection for the fast deblurring and
% interpolation method. This function implements the gradient of the
% level one norm between the blurred version of the current deblurred HR estimate and
% the original blurred HR estimate created by the median and shift method.
%
% Inputs:
% Xn - The current estimate of the deblurred HR image
% Z -  The original blurred estimate of the HR image
% A -  Normilzation factor for each pixel
% Hpsf - The PSF function (common to all frames and space invariant)
%
% Outpus:
% The backprojection of the sign of the residual error
function G=FastGradientBackProject(Xn, Z, A, Hpsf) 

% Blur the current HR estimate
Zn = imfilter(Xn, Hpsf, 'symmetric');

% Deblur the normalized sign of Xdiff
Gsign = sign(A.*(Zn-Z));

% Unblur the backprojected image
G = imfilter(A.*Gsign, flipud(fliplr(Hpsf)), 'symmetric');
