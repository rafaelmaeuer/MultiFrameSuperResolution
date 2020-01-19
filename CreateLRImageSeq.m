% This is a utility function which takes an HR image and generates a LR
% image sequence after motion, blurring and decimation and noise addition
% operations are performed.
%
% Inputs:
% HR - The original HR image
% resFactor - The resolution increment factor
% noiseSig - The sigma of the white gaussian noise to be added to the LR images.
% Hpsf - The PSF function (common to all frames and space invariant)
% N - The number of frames in the LR sequence to generate
%
% Outputs:
% LR - The sequence of LR images
% D  - The motion for each frame
function [LR, D]=CreateLRImageSeq(HR, resFactor, noiseSig, Hpsf, N)

% Blur the HR image
Z = imfilter(HR, Hpsf, 'symmetric');

% Randomize shifts
D=1+floor(rand(N,2)*(resFactor*2-1));

LR = zeros([floor(size(HR)/resFactor)-2 N]);

% Sample Z based on shift
for i=1:N

  LR(:,:,i)=Z(D(i,2):resFactor:D(i,2)+size(LR,1)*resFactor-1, D(i,1):resFactor:D(i,1)+size(LR,2)*resFactor-1);
  
end

% Add white guassian noise.
LR = LR + randn(size(LR))*noiseSig;