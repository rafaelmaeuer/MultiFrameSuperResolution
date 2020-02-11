% Image-Registration of a given sequence of images
%
% This function uses a simplified Lucas-Kanade method to
% register images solely based on their translatory
% deviation. It uses a hierarchical gradient-based
% optimization method, using 6 levels of Low-Pass filtering
% for calculation of the optical flow parameters between
% two images
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
function [LR_reg, Tvec, iter, err]=RegisterImageSeq(app, stack)

    % Init variables
    iter = 0; err = 0;
    
    % Get baseframe
    baseFrame = squeeze(stack(:,:,1));
    height = size(baseFrame,1);
    width = size(baseFrame,2);
    
    % Create the region of continuous flow for lowest hierarchy
    % level of the gaussian pyramid
    roi=[2 2 size(stack,1)-1 size(stack,2)-1];

    for i=2:size(stack,3)
        ShowProgress(app, 'Image Registration in progress...', (i/size(stack,3))*100);

        % Register current image to previous frame
        [D, k, e] = PyramidalLKOpticalFlow(baseFrame, squeeze(stack(:,:,i)), roi);

        % Project the 2D-motion vector onto the general affine
        % transformation matrix
        Tvec(i,:) = D;

        % Perform Image registration
        tform = affine2d([ 1 0 0; 0 1 0; D(1) D(2) 1]);
        I = imwarp(stack(:,:,i),tform,'cubic','FillValues',128);
        LR_reg(:,:,i) = APP_cropImage(I, width, height);

        % Sum up the iterations and errors
        iter = iter + k;
        err = err + e;
    end
  
end