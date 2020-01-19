% Computes the gradient of the regularization term of the super-resolution
% minimization function. This function implements the gradient of the
% bilateral-filter.
%
% Inputs:
% Xn - The current estimate of the denoise/deblured image
% P  - The spatial window size (radius)
% alpha - The exponential decay coefficient
%
% Outpus:
% The backprojection of the sign of the residual error
function G=GradientRegulization(Xn, P, alpha)

G=zeros(size(Xn));

% Create an inflated version of Xn so shifting operation is simpler
Xpad = padarray(Xn, [P P], 'symmetric');

% Compute a grid of l=-P:P and m=0:P such that l+m>=0
for l=-P:P
  for m=-P:P

    % Shift HR by l and m
    Xshift = Xpad(1+P-l:end-P-l, 1+P-m:end-P-m);

    % Subtract from HR image and compute sign
    Xsign = sign(Xn-Xshift);
    
    % Shift Xsign back by -l and -m
    Xsignpad = padarray(Xsign, [P P], 0);

    Xshift = Xsignpad(1+P+l:end-P+l, 1+P+m:end-P+m);

    G = G + alpha.^(abs(l)+abs(m)).*(Xsign-Xshift);

  end
end