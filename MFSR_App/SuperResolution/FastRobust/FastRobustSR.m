% Implements the fast and robust super-resolution method. This funtion
% first compute an estimation of the blurred HR image, using the median and
% shift method. It then uses the bilateral filter as a regulating term
% for the debluring and interpolation step.
%
% Inputs:
% app - instance of the main app
% LR - A sequence of low resolution images
% Tvec - The tranlational motion for each LR frame
% resFactor - The resolution increment factor
% Hpsf - The PSF function (common to all frames and space invariant)
% params - property structure used to control the algorithm parameters
%
% Outputs:
% HR - The estimated HR image
% iter - Steps needed for estimation
function [HR, iter]=FastRobustSR(app, LR, Tvec, resFactor, Hpsf, params)

    % project the translation to the new image size, rounded to
    % the nearest neighbour
    D = round(Tvec.*resFactor);

    % backproject the rounded vector to the initial size
    Dr = floor(D/resFactor);
    D = mod(D,resFactor)+resFactor;

    [X,Y] = meshgrid(1:size(LR, 2), 1:size(LR, 1));

    for i=1:size(LR, 3)
        LR(:,:,i) = interp2(X+Dr(i,1), Y+Dr(i,2), LR(:,:,i), X, Y, '*nearest');
    end
    stack_r = LR(3:end-2,3:end-2,:);

    % Compute initial estimate of blurred HR by the means of MedianAndShift
    [Z, A] = MedianAndShift(stack_r, D, [(size(stack_r,1)+1)*resFactor-1 (size(stack_r,2)+1)*resFactor-1], resFactor);

    % Deblur the HR image and regulate using bilatural filter

    % Loop and improve HR in steepest descent direction
    HR = Z;
    iter = 1;

    while iter<params.maxIter
        ShowProgress(app, ' Estimating High Resolution image', (iter/params.maxIter*100));

        % Compute gradient of the energy part of the cost function
        Gback = FastGradientBackProject(HR, Z, A, Hpsf);

        % Compute the gradient of the bilateral filter part of the cost function
        Greg = GradientRegulization(HR, params.P, params.alpha);

        % Perform a single SD step
        HR = HR - params.beta.*(Gback + params.lambda.* Greg);

        iter=iter+1; 
    end

end
