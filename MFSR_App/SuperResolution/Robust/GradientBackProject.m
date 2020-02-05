% Computes the gradient backprojection for the super-resolution
% minimization function. This function implements the gradient of the
% level one norm between the projection of the estimated HR image and each
% of the LR images.
%
% Inputs:
% Xn - The current estimate of the HR image
% LR - A sequence of low resolution images
% Tvec - The tranlational motion for each LR frame
% Hpsf - The PSF function (common to all frames and space invariant)
% Dres - The resolution increment factor
%
% Outpus:
% The backprojection of the sign of the residual error
% TODO: rename Tvec to D (project uniformity?)
function G=GradientBackProject(Xn, LR, Tvec, Hpsf, Dres) 

% Note that shift and blur are comutative, so to improve runtime, we first
% filter the HR image
Zn = imfilter(Xn, Hpsf, 'symmetric');

% Allocate shifted and decimated HR image
HRsd = zeros(size(LR));

for k=1:size(LR,3)

  % Shift and decimate HR image for each frame k
  % TODO: find meaningful variable names
%   temp1=(Tvec(k,2):Dres:(size(LR,1)-1))*Dres; % (dY:H_image) in steps of resFactor
%   temp2= Tvec(k,2); % resFactor + dY
%   temp3=(Tvec(k,1):Dres:(size(LR,2)-1))*Dres; % (dX:W_image) in steps of resFactor
%   temp4= Tvec(k,1); % resFactor + dx
%   temp5= temp1+temp2;
%   temp6= temp3+temp4;
%   HRsd(:,:,k)=Zn(temp5, temp6);
  % TODO: Find cause for Error when running Robust Algorithm first
  HRsd(:,:,k)=Zn(Tvec(k,2):Dres:(size(LR,1)-1)*Dres+Tvec(k,2),Tvec(k,1):Dres:(size(LR,2)-1)*Dres+Tvec(k,1));
  
end

% Compute the sign between HRsd-LR
Gsign = sign(HRsd-LR);

HRsd = zeros([size(Xn) size(LR,3)]);

% Back project Gsign to HR space
for k=1:size(LR,3)

  % Upsample and shift LR sign image for each frame k
  % TODO: split variables and give it meaningful names?
  HRsd(Tvec(k,2):Dres:(size(LR,1)-1)*Dres+Tvec(k,2),Tvec(k,1):Dres:(size(LR,2)-1)*Dres+Tvec(k,1),k)=Gsign(:,:,k);
  
end

% Unblur the backprojected image
G = imfilter(HRsd, flipud(fliplr(Hpsf)), 'symmetric');

% Compute the sum over k of the backprojected gradient
G = sum(G, 3);
