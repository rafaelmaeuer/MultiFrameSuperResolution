% Implements a simple cubic-spline interpolation of a single image. This
% image is then deblured using the same method as in the Fast and Robust
% method.
%
% Inputs:
% app - instance of the main app
% LR - A sequence of low resolution images
% resFactor - The resolution increment factor
% Hpsf - The PSF function (common to all frames and space invariant)
% params - property structure used to control the algorithm parameters
%
% Outputs:
% HR - The estimated HR image
% iter - Steps needed for estimation
function [HR, iter]=CubicSplineInterp(app, LR, resFactor, Hpsf, params)

    % Initialize guess as interpolated version of LR
    [X,Y]=meshgrid(0:resFactor:(size(LR,2)-1)*resFactor, 0:resFactor:(size(LR,1)-1)*resFactor);
    [XI,YI]=meshgrid(resFactor+1:(size(LR,2)-2)*resFactor-1, resFactor+1:(size(LR,1)-2)*resFactor-1);

    Z=interp2(X, Y, squeeze(LR(:,:,1)), XI, YI, '*spline');

    % Deblur the HR image and regulate using bilatural filter

    % Loop and improve HR in steepest descent direction
    HR = Z;
    iter = 1;
    A = ones(size(HR));

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
