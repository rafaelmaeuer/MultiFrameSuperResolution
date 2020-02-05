function [Z, D, A] = APP_prepareRSR(LR, resFactor, Tmat)
    % Preparation
    stack = LR(3:end-2,3:end-2,:);
    % Extract the motion vectors from the transformation
    % matrices and round translation to nearest neighbor
    % makes the translation pixel-perfect
    Dr = round(squeeze(Tmat(1:2, 3,:)).*resFactor).';
    D = mod(Dr,resFactor)+resFactor;

    % Compute initial estimate of blurred HR by the means of MedianAndShift
    [Z, A]=MedianAndShift(stack, D, [(size(stack,1)+1)*resFactor-1 (size(stack,2)+1)*resFactor-1], resFactor);

end

