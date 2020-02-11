% Matlab Image Registration
%
% Inputs:
% app - instance of the main app
% stack - A sequence of low resolution images
%
% Outputs:
% LR_reg - The low resolution stack with registered images
% Tvec - The tranlational motion for each LR frame
function [LR_reg, Tvec]=RegisterImageSeqMatlab(app, stack)
    
    % Get baseframe
    baseFrame = squeeze(stack(:,:,1));
    height = size(baseFrame,1);
    width = size(baseFrame,2);
    
    % Init matlab image registration
    [optimizer, metric] = imregconfig('monomodal');

    % Iterate through all frames except the base frame
    for i=2:size(stack, 3)
        ShowProgress(app, 'Image Registration in progress...', (i/size(stack,3)*100));

        % Get transformation matrix from matlab image registration
        tform = imregtform(stack(:,:,i), baseFrame, 'affine', optimizer, metric);

        % Extract the translation vector and set translation zero
        Tvec(i,:) = tform.T(3,1:2);
        tform.T(3,1:2) = [0 0];

        % Warp the current frame using the calculated transformation matrix
        I = imwarp(stack(:,:,i), tform, 'cubic', 'FillValues', 128);

        % Crop Image to initial size
        LR_reg(:,:,i) = CropImage(I, width, height);
    end
    
end

