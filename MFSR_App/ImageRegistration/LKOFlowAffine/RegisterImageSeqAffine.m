% Affine Image-Registration of a given sequence of images
%
% This function uses a Lucas-Kanade method to register images 
% based on their translatory, rotational and shearing deviation.
% It uses a hierarchical gradient-based optimization method, 
% using 6 levels of Low-Pass filtering for calculation of the 
% optical flow parameters between two images
%
% Inputs:
% app - instance of the main app
% stack - A sequence of low resolution images
%
% Outputs:
% LR_reg - The low resolution stack with registered images
% Tvec - The tranlational motion for each LR frame
% iter - Steps needed for registration
% err - Sum of error during registration
function [LR_reg, Tvec, iter, err]=RegisterImageSeqAffine(app, stack)

    % Init variables
    iter = 0; err = 0;
    
    % Get baseframe
    baseFrame = squeeze(stack(:,:,1));
    height = size(baseFrame,1);
    width = size(baseFrame,2);
    
    % Create the region of continuous flow for lowest hierarchy
    % level of the gaussian pyramid
    roi=[2 2 size(stack,1)-1 size(stack,2)-1];

    % Initialize transformation matrix
    D = [1 0 0; 0 1 0];

    for i=2:size(stack,3)

        ShowProgress(app, 'Image Registration in progress...', (i/size(stack,3))*100);

        % Register current image to previous frame
        dc = PyramidalLKOpticalFlowAffine(baseFrame, squeeze(stack(:,:,i)), roi);

        % Set the current frame as base-frame for the next iteration
        baseFrame = squeeze(stack(:,:,i));

        % Add current displacement dc to D (This is actually concatinating the two
        % affine matrixes)
        D = D + reshape(dc, 2, 3)*(eye(3)+[D;0 0 0]);

        % Compute displacement at current level
        [D,k,e] = IterativeLKOpticalFlowAffine(squeeze(stack(:,:,1)), squeeze(stack(:,:,i)), roi, D);

        % Set the return value for translation vector
        Tvec(i,:) = [D(1,3), D(2,3)];

        % Perform MATLAB Image Registration
        tform = affine2d([ D(1,1) D(2,1) 0; D(1,2) D(2,2) 0; 0 0 1]);
        I = imwarp(stack(:,:,i),tform,'cubic','FillValues',128);

        % LKFlowAffine performs a weird zoom, so we have to
        % resize the images to get size-normalized results
        LR_reg(:,:,i) = imresize(I, [height, width]);

        % Sum up the iterations and errors
        iter = iter + k;
        err = err + e;
    end

end