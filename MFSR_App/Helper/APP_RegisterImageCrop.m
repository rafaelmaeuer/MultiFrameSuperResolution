function [imOut] = APP_RegisterImageCrop(theImage,theTransform)

% Compute coordinates for resampling
[X, Y] = meshgrid(1:size(theImage,2), 1:size(theImage, 1));

% calculate projection to new position
D = theTransform *[X(:)';Y(:)' ; ones(1,length(X(:)))];

% Update coorindates with displacement
X = reshape(D(1,:), size(X));
Y = reshape(D(2,:), size(Y));

% calculate the interpolated version of the new image
imOut = interp2(theImage,X,Y);

figure; image(imOut)

end

