% Pyramid = GuassianPyramid(img, levels)
%
% Computes a gaussian pyramid of the given number levels.
% The image is downsampled by 2 in each level. Before downsampling an
% antialiasing filter using a kernel: [1 4 6 4 1]/16 is performed.
%
% Inputs:
% img - The base image of the highest level in the pyramid
% levels - The number of levels to be used. The parameter is optional.
%                     A default level of 4 will be used of parameter is
%                     ignored.
% Outpus:
% A cell array  of images. Index 1 will include the input img. Index
% k=levels will include the lowest resolution image.
function pyramid = GuassianPyramid(img, levels)

if nargin<2
  levels=4;
end

% Initialize pyramid
pyramid{1}=img;

% Build pyramid
for i=2:levels
  pyramid{i}=GaussianDownSample(pyramid{i-1});
end