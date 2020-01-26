% Computes the gradient backprojection for the super-resolution
% minimization function. This function implements the gradient of the
% level one norm between the projection of the estimated HR image and each
% of the LR images.
%
% Inputs:
% Xn - The current estimate of the HR image
% LR - A sequence of low resolution images
% Fmot - The tranlational motion for each LR frame
% Hpsf - The PSF function (common to all frames and space invariant)
% Dres - The resolution increment factor
%
% Outpus:
% The backprojection of the sign of the residual error
function G=GradientBackProject(Xn, LR, Fmot, Hpsf, Dres) 

% Note that shift and blur are comutative, so to improve runtime, we first
% filter the HR image
Zn = imfilter(Xn, Hpsf, 'symmetric');

% Allocate shifted and decimated HR image
HRsd = zeros(size(LR));

for k=1:size(LR,3)

  % Shift and decimate HR image for each frame k
  HRsd(:,:,k)=Zn(Fmot(k,2):Dres:(size(LR,1)-1)*Dres+Fmot(k,2),Fmot(k,1):Dres:(size(LR,2)-1)*Dres+Fmot(k,1));
  
end

% Compute the sign between HRsd-LR
Gsign = sign(HRsd-LR);

HRsd = zeros([size(Xn) size(LR,3)]);

% Back project Gsign to HR space
for k=1:size(LR,3)

  % Upsample and shift LR sign image for each frame k
  HRsd(Fmot(k,2):Dres:(size(LR,1)-1)*Dres+Fmot(k,2),Fmot(k,1):Dres:(size(LR,2)-1)*Dres+Fmot(k,1),k)=Gsign(:,:,k);
  
end

% Unblur the backprojected image
G = imfilter(HRsd, flipud(fliplr(Hpsf)), 'symmetric');

% Compute the sum over k of the backprojected gradient
G = sum(G, 3);
