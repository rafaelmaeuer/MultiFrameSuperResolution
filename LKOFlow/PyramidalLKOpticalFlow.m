% d=PyramidalLKOpticalFlow(img1, img2, roi)
%
% Computes the Lucas-Kanade optical flow using a gaussian pyramid.
% The function builds a pyramid for images 1 and 2 and computes the flow
% in each level, and using the flow an lower levels as initial estimates
% for the flow in higher levels. The flow is computed as the displacement
% from img2 to img1.
%
% Inputs:
% img1 - The source image.
% img2 - The destination image.
% roi  - A region in img1 in which the flow should be computed.
%
% Outpus:
% d - The computed displacement
function  d=PyramidalLKOpticalFlow(img1, img2, roi)

img1=double(img1);
img2=double(img2);

roiSize = 1+roi(3:4)-roi(1:2);

% Compute boundry for number of levels. We define two bounds. K<=6. In most
% cases more than 6 levels is not needed => the image becomes too small.
% Also, the roi at level k should not be too small. We define here that roi
% should not become less than 4x4 pixels.
%
% Since the size of the ROI at level k is Roiwidth/2^k we have
% Roiwidth/2^k >= 4 => 2^k<=Roiwidth/4 => k<log2(Roiwidth/4) =
% log2(Roiwidth)-2

levels = min([floor(log2(min(roiSize))-2), 6]);

% Compute Pyramids for image 1 and 2
pyramid1 = GuassianPyramid(img1, levels);
pyramid2 = GuassianPyramid(img2, levels);

d=[0 ; 0];

% Start with initial translation 0 at lowest pyramid
for l=levels:-1:1

  % Transform displacement to current level
  d=d*2;

  % Compute current location of ROI
  tl=ceil(roi(1:2)./2^(l-1));

  % Compute current size of ROI
  sz=floor(roiSize./2^(l-1));

  % Make sure origin is at least 2
  tl=max(tl, [2 2]);
  br=min(tl+sz-1, size(pyramid1{l})-1);
  
  % Compute displacement at current level
  d=IterativeLKOpticalFlow(pyramid1{l}, pyramid2{l}, [tl br], d);

end




