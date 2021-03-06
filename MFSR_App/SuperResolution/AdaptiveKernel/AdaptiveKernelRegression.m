% This program implements the algorithm presented in the following paper
% Mohammad Moinul Islam, Vijayan K. Asari, Mohammed Nazrul Islam, and 
% Mohammad A. Karim, "Video super-resolution by adaptive kernel regression" 
% Lecture Notes in Computer Science, Published by Springer-Verlag Berlin/
% Heidelberg ISSN: 0302-9743 (Print) 1611-3349 (Online), Advances in Visual
% Computing, Edited by G. Bebis et al. Proceedings of the International 
% Symposium on Visual Computing -ISVC - 2009: ISBN: 978-3-642-10519-7, 
% vol. 5876/2009, pp. 799-806, November 2009
%
% Inputs:
% app - instance of the main app
% LR - A sequence of low resolution images
% resFactor - The resolution increment factor
% params - property structure used to control the algorithm parameters
%
% Outputs:
% HR - The estimated HR image
% iter - Steps needed for estimation
function [HR, iter]=AdaptiveKernelRegression(app, LR, resFactor, params)

% set input buffer
buffer = LR;

% get input dimensions
[m,n,s] = size(LR);

% what are these for?
h= 0.5;
hr= 255;

% set scaling factor
f= resFactor;
% set position by scaling factor
fp = f-1;
% set iterator by scaling factor
fi= f^2;

% waitbar
iter= 1;
maxIter= m*n;

% loop through y (image height) minus scaling factor as border
for i= 1:m-(f)
    ShowProgress(app, ' Estimating High Resolution image', (iter/maxIter*100));

    % loop through x (image width) minus scaling factor as border
    for j= 1:n-(f)

        % init arrays for kernel-regression and kernel-regression-frame (s) pixel values
        for l= 1:fi
            kr(l) = 0;
            krf(l)= 0;
        end

        % loop through s frames or maximum number of iterations
        if(params.maxIter < s)
            s = params.maxIter;
        end

        for k= 1:s

            % create Gaussian kernel
            kr_bilat= (1/(f*pi*hr*hr))*exp(-(norm(buffer(i,j,k)-mean(buffer(i,j,:))).^2)/(f*hr*hr));

            % init line and index values
            line = 0; idx = 0;
            for o= 1:fi

                % sum kernel-regression and kernel-regression-frame (k) pixel values
                kr(o)= kr(o) + (1/(f*pi*h*h))*exp(-(norm((i+line-i)^2+(j+idx-j)^2).^2)/(f*h*h))*kr_bilat;
                krf(o)= krf(o) + buffer(i+line,j+idx,k)*(1/(f*pi*h*h))*exp(-(norm((i+line-i)^2+(j+idx-j)^2).^2)/(f*h*h))*kr_bilat;
                idx = idx + 1;

                % if index is multiple of resolution factor, increment line
                if(mod(o, f) == 0)
                    line = line + 1;
                    idx = 0;
                end
            end

        end

        % init line and index values
        line = 1; idx = 1;
        for r= 1:fi

            % fill temporary quadratic matrix in size of resolution factor
            tmp(line, idx)= krf(r)/kr(r);
            idx = idx + 1;

            % if index is multiple of resolution factor, increment line
            if(mod(r, f) == 0)
                line = line + 1;
                idx = 1;
            end
        end

        % write pixel values from temporary matrix to correct position
        ynew(f*i-fp:f*i,f*j-fp:f*j)= tmp;
        iter=iter+1;
    end
    iter=iter+1;
end

% median filter
HR= medfilt2(ynew,[f f]);
