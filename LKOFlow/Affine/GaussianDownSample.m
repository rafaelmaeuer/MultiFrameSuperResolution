% img=GaussianDownSample(img, kernel)
%
% Downsamples the given image by 2. A lowpass filter is performed prior to 
% the downsampling.
%
% Inputs:
% img - The image to downsample
% kernel - An optional kernel for downsampleing. If the kernel is not
% provided, the default seperable filter [1 4 6 4 1]/16 is used.
%
% Outputs:
% The downsampled image
function img=GaussianDownSample(img, kernel)

if nargin<2
    kernel = [1 4 6 4 1]/16;
end

if size(kernel, 1) ==1||size(kernel, 2) ==1
    seperable=1;
end

% Perform antialiasing filter
img=imfilter(img, kernel, 'symmetric');

% If filter is seperable perform traposed filter too
if seperable
  img=imfilter(img, kernel', 'symmetric');
end

% Down sample
img = img(1:2:end, 1:2:end);