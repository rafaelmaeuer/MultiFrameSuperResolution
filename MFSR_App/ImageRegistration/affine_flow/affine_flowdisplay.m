function affine_flowdisplay(flow, im, step, col)
%AFFINE_FLOWDISPLAY Display affine flow field.
%   AFFINE_FLOWDISPLAY(FLOW, IM) takes FLOW, a structure returned by
%   AFFINE_FLOW and an image IM. Vectors are drawn on the image
%   representing the flow at selection of points. The point at which the
%   flow is computed is marked with a circle and a flow vector is shown as
%   a line pointing away from the circle.
%
%   AFFINE_FLOWDISPLAY(FLOW, REGION) where REGION is a row vector, draws
%   the vectors in the current figure without setting a background. REGION
%   can have 4 elements and specifies the region in which to draw vectors
%   in the format [Ymin Ymax Xmin Xmax] (for consistency with AFFINE_FLOW).
%   REGION may have 2 elements giving [Ymax Xmax]; Ymin and Xmin both
%   default to 1. Thus AFFINE_FLOWDISPLAY(FLOW, SIZE(IM)) displays the
%   flow for an image on a blank canvas.
%
%   AFFINE_FLOWDISPLAY(..., STEP) steps the separation between vector
%   locations. If STEP is a scalar it sets the separation on both axes;
%   otherwise it must be a vector with the form [STEPY STEPX], or an empty
%   matrix to use the default setting of about 10 vectors across the image.
%
%   AFFINE_FLOWDISPLAY(..., COL) sets the colour in which to draw the
%   vectors. COL must be a single colour character, as for PLOT.
%
%   See also: AFFINE_FLOW, AFFINE_FLOWDEMO

% Copyright David Young 2010

if isequal(size(im), [1 2])
    im = [1 im(1) 1 im(2)];
end
if isequal(size(im), [1 4])
    r0 = im([3 1]);
    r1 = im([4 2]);
    axis(im([3 4 1 2]));
    axis equal;
    set(gca,'YDir','reverse');
else
    r0 = [1 1];
    r1 = [size(im,2) size(im,1)];
    imshow(im, []);
end

if nargin < 3 || isempty(step)
    step = (r1 - r0)/10;
elseif isscalar(step)
    step = [step step];
else
    step = [step(2) step(1)];   % for consistency with region
end
    
if nargin < 4
    col = 'g';
end

hold on;
disp_vecs(r0, step, r1, flow, col);
hold off;

end

%-----------------------------------------------------------------------

function [xs1, ys1, xs2, ys2] = flow2offsets(reg0, step, reg1, flow)
% Generate a set of position changes for display on a rectangular grid

% get starting positions
size = reg1 - reg0;
nsteps = floor(size./step);
starts = reg0 + (size - step.*nsteps)/2;
xs1 = linspace(starts(1), starts(1)+nsteps(1)*step(1), nsteps(1)+1);
ys1 = linspace(starts(2), starts(2)+nsteps(2)*step(2), nsteps(2)+1);
% expand to 2D grid
[xs1, ys1] = meshgrid(xs1, ys1);
xs1 = xs1(:);
ys1 = ys1(:);

% compute end positions by applying flow
w = affine_flow.warp(flow);
pos2 = [xs1 ys1 ones(length(xs1),1)] * w;
xs2 = pos2(:,1);
ys2 = pos2(:,2);

end

%-----------------------------------------------------------------------

function disp_vecs(reg0, step, reg1, flow, colour)
% Display the flow field as a set of vectors on a rectangular grid

[xs1, ys1, xs2, ys2] = flow2offsets(reg0, step, reg1, flow);

plot(xs1, ys1, [colour 'o']);
plot([xs1 xs2].', [ys1 ys2].', [colour '-']);

end
