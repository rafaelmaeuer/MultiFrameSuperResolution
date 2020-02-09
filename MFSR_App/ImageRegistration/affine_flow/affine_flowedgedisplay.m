function affine_flowedgedisplay(flow, im1, im2)
%AFFINE_FLOWEDGEDISPLAY displays image match under affine flow
%   AFFINE_FLOWEDGEDISPLAY(FLOW, IM1, IM2) takes FLOW, a structure returned
%   by AFFINE_FLOW, and the two images that were used to estimate it. The
%   edges from the images are displayed in the current figure, using the
%   following colours:
%
%       green: edges from IM1
%       blue: edges from IM2
%       red: edges from IM1 after warping by the flow field
%
%   Good results are indicated if the red edges are close to the blue
%   edges.
%
%   See also: AFFINE_FLOW, AFFINE_FLOWDEMO

% Copyright David Young 2010

w = affine_flow.warp(flow);
t = maketform('affine', w);
im1trans = imtransform_same(im1, t);

e1 = edge(im1, 'canny');
e1trans = edge(im1trans, 'canny');
e2 = edge(im2, 'canny');

% combine the edges to show them in different colours
edges(:,:,1) = 1 - (e1 | e2);
edges(:,:,2) = 1 - (e2 | e1trans);
edges(:,:,3) = 1 - (e1 | e1trans);

imshow(edges);

end

