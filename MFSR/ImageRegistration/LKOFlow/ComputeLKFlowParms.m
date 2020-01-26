% [Ht, G]=ComputeLKFlowParms(img)
%
% Computes the optical flow parameters. The image is derived in DX and DY
% directions and matrix G is computed.
%
% Inputs:
% img - The image to compute optical flow parameters. Note that the image
%       must contain a border of size 1 since the DX and DY operation work
%       on 3x3 regions.
%
% Outpus:
% Ht - The transpose of matrix H, which defines the LK affine approximation of an shifted
% image
%
% G  - The matrix G=H'*H
%
function [Ht, G]=ComputeLKFlowParms(img)

% Compute Ix and Iy derivatives using sobel operator
Hdy = fspecial('sobel');
Hdx = Hdy';

Ix=imfilter(img, Hdx);
Iy=imfilter(img, Hdy);

Ix=Ix(2:end-1,2:end-1);
Iy=Iy(2:end-1,2:end-1);

% Compute G
H = [Ix(:) Iy(:)];
Ht=H';
G=Ht*H;
