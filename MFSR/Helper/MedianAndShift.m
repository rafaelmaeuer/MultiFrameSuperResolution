% [Z,A]=MedianAndShift(LR, D, HRsize, Dres)
%
% Estimates (Robustly) the blurred high resolution as a median of the LR
% images after upsampleing and shifting to correct position. This is
% prooved by Farisu et al in "Fast and Robust Multiframe Super Resolution"
% to be the estimator (with no regularization) which minimizes the L1
% norm. 
%
% An efficient implementation is performed here in which we note that the
% high resolution image consists of the median of all the LR images which
% have the same displacement value. That is, the LR images are partitioned
% to r^2 unique groups, each group have equal displacement values in the HR
% image. The median of each group is then upsampled and shifted into the HR
% image. 
% 
% In some cases, not all r^2 displacements exists. This leaves us with
% undetermined cases. In these case some other interpolatoin method needs
% to be implemented. In this implementation we "fill" the holes using
% a spatial median filter.
%
% Inputs:
%
% LR - The sequence of low resolution images
% D - The displacement vector for each frame
% HRsize - The size of the HR image
% Dres - The resolution scale factor.
%
% Outputs:
%
% Z - The original blurred estimate of the HR image
% A - Normilzation factor for each pixel

function [Z,A]=MedianAndShift(LR, D, HRsize, Dres)

% Allocate high resolution image
Z = zeros(HRsize);
A = ones(HRsize);

S = zeros(Dres);

% Loop on each possible displacement value (should be much less than the
% number of images in LR for over determined solution)
for x=Dres:2*Dres-1
  for y=Dres:2*Dres-1
    
    % Find displacement values
    % I is a vector of length(number of frames) which contains boolean
    % values. An entry of I is true, if the displacement matches the
    % iteration indices.
    I = D(:,1)==x & D(:,2)==y;
    
    % returns the number of valid displacements
    len = length(find(I==true));
    
    % Handle only cases in which there is at least one LR at this
    % displacement
    if len>0

      % Indicate data exists for this shift by writing '1' at the
      % corresponding pixel location in S
      S(x-Dres+1, y-Dres+1)=1;

      % Fill Matrices Z and A by calculating median
      Z(y:Dres:size(Z, 1),x:Dres:size(Z, 2))= median(LR(:,:,I), 3);
      A(y:Dres:size(Z, 1),x:Dres:size(Z, 2))=len;
    end

  end
end

% Find under-determined shifts
[X,Y] = find(S==0);
% 
if ~isempty(X)
  
  % Compute a median filter with window size of the resolution factor
  % assuming that more than 50% of the shifts should be determined 
  Zmed=medfilt2(Z, [Dres Dres]);
  
  % Loop on each hole and fill with median
  for i=1:length(X)
    x =X(i)+Dres-1;
    y =Y(i)+Dres-1;

    Z(y:Dres:size(Z, 1),x:Dres:size(Z, 2))=Zmed(y:Dres:size(Z, 1),x:Dres:size(Z, 2));
    
  end
  
end

A = sqrt(A);
