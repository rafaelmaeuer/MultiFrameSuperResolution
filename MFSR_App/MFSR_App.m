classdef MFSR_App < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        MFSRToolUIFigure          matlab.ui.Figure
        IMG_LR                    matlab.ui.control.Image
        IMG_HR                    matlab.ui.control.Image
        LABEL_heading_LRstack     matlab.ui.control.Label
        LABEL_heading_HRstack     matlab.ui.control.Label
        PANEL_imInfo              matlab.ui.container.Panel
        GRID_imInfo               matlab.ui.container.GridLayout
        LABEL_imInfo_nframes      matlab.ui.control.Label
        LABEL_imInfo_imHeight     matlab.ui.control.Label
        LABEL_imInfo_format       matlab.ui.control.Label
        LABEL_imInfo_imWidth      matlab.ui.control.Label
        VAL_imInfo_nframes        matlab.ui.control.Label
        VAL_imInfo_format         matlab.ui.control.Label
        VAL_imInfo_imWidth        matlab.ui.control.Label
        VAL_imInfo_imHeight       matlab.ui.control.Label
        BTN_next                  matlab.ui.control.Button
        BTN_prev                  matlab.ui.control.Button
        BTN_loadFile              matlab.ui.control.Button
        PANEL_bench               matlab.ui.container.Panel
        GRID_bench                matlab.ui.container.GridLayout
        LABEL_heading_IR          matlab.ui.control.Label
        LABEL_heading_SR          matlab.ui.control.Label
        GRID_IR_bench             matlab.ui.container.GridLayout
        LABEL_IR_time             matlab.ui.control.Label
        LABEL_IR_iter             matlab.ui.control.Label
        LABEL_IR_err              matlab.ui.control.Label
        VAL_IR_t                  matlab.ui.control.Label
        VAL_IR_n                  matlab.ui.control.Label
        VAL_IR_err                matlab.ui.control.Label
        GRID_SR_bench             matlab.ui.container.GridLayout
        LABEL_SR_time             matlab.ui.control.Label
        LABEL_SR_iter             matlab.ui.control.Label
        LABEL_SR_err              matlab.ui.control.Label
        VAL_SR_t                  matlab.ui.control.Label
        VAL_SR_n                  matlab.ui.control.Label
        VAL_SR_err                matlab.ui.control.Label
        BTN_enhance               matlab.ui.control.Button
        RADIOGROUP_IR_Method      matlab.ui.container.ButtonGroup
        RADIO_IR_MATLAB           matlab.ui.control.RadioButton
        RADIO_IR_LKFlowAffine     matlab.ui.control.RadioButton
        RADIO_IR_LKFlow           matlab.ui.control.RadioButton
        RADIOGROUP_SR_Method      matlab.ui.container.ButtonGroup
        RADIO_SR_KernelReg        matlab.ui.control.RadioButton
        RADIO_SR_Cubic            matlab.ui.control.RadioButton
        RADIO_SR_Robust           matlab.ui.control.RadioButton
        RADIO_SR_FastRobust       matlab.ui.control.RadioButton
        BTN_reset                 matlab.ui.control.Button
        PANEL_SR_params           matlab.ui.container.Panel
        ResfactorLabel            matlab.ui.control.Label
        PARAM_ResFactor           matlab.ui.control.NumericEditField
        IterationsEditFieldLabel  matlab.ui.control.Label
        PARAM_Iterations          matlab.ui.control.NumericEditField
        BTN_saveLRFrame           matlab.ui.control.Button
        BTN_saveHRFrame           matlab.ui.control.Button
        StatusPanel               matlab.ui.container.Panel
        GridLayout                matlab.ui.container.GridLayout
        perc0                     matlab.ui.control.Label
        perc1                     matlab.ui.control.Label
        perc2                     matlab.ui.control.Label
        perc3                     matlab.ui.control.Label
        perc4                     matlab.ui.control.Label
        perc5                     matlab.ui.control.Label
        perc6                     matlab.ui.control.Label
        perc7                     matlab.ui.control.Label
        perc8                     matlab.ui.control.Label
        perc9                     matlab.ui.control.Label
    end

    
    properties (Access = private)
        LRstack     % Container for the low resolution video frames
        LR_reg
        TM          % the affine transformation matrix
        imReg_flag  % Flag that indicates if image registration was already done once
        HRimage     % HR image
        frameCnt    % counter to track the currently displayed frame
    end
    
    methods (Access = private)
        
        function [params, resFactor, Hpsf] = collectParams(app)
            
            % Parameter for resolution factor
            resFactor = app.PARAM_ResFactor.Value;

            % Parameter for gaussian filter
            psfSize = 3;
            psfSig = 1;
            Hpsf = fspecial('gaussian', [psfSize psfSize], psfSig);

            % Parameter for image registration
            params.alpha = 0.7;     % Spatial correction Factor for regularization term
                                    % of the SR optimization
            params.beta = 1;        % Step Size of the steepest descent optimization
            params.lambda = 0.04;   % Regularization Parameter weighting the 
                                    % regularization cost-function
                        
            params.P = 4;           % Number of shifts to calculate the regularization
                                    % cost-function
                        
            % stop criterion (max numbers of iterations in  steepest descent optimization)
            params.maxIter = app.PARAM_Iterations.Value;
        end
        
        function [Tmat, iter, err] = getTransform(app)
            % Returns a 3x3xF stack of 3x3 image transformation matrices for each
            % frame f with respect to the reference frame (first frame in
            % stack).
            showprogress(app, 'Image Registration in progress...', 0);
            stack = app.LRstack;
            Tmat = zeros(3,3,size(stack,3));
            iter = 0;
            err = 0;
            baseFrame = squeeze(stack(:,:,1));
            
            Tmat(:,:,1) = eye(3,3);
            
            % check if MATLAB registration Mode was selected
            if app.RADIO_IR_MATLAB.Value
                
                % Init matlab image registration
                [optimizer, metric] = imregconfig('monomodal');
            
                % Iterate through all frames
                for i=2:size(stack, 3)
                    showprogress(app, 'Image Registration in progress...', (i/size(stack,3)*100));
                    % Get transformation vector from matlab image registration
                    Tmat_i = imregtform(app.LRstack(:,:,i), baseFrame, 'affine', optimizer, metric);
                    
                    % For some reason, the translation parameters are set
                    % in the wrong position by the matlab method, so we
                    % have to make some arrangements
                    Tmat(:,:,i) = Tmat_i.T.';
                    Tmat(1:2,1:2,i) = Tmat_i.T(1:2,1:2);
                end
            end
            
            if app.RADIO_IR_LKFlowAffine.Value
                
                roi=[2 2 size(stack,1)-1 size(stack,2)-1];
                
                Mprev = squeeze(stack(:,:,1));
                
                % Tmat(:,:,1) = eye(3,3);
                D = Tmat(1:2,1:3,1);
                
                for i=1:size(stack,3)
                    
                    showprogress(app, 'Image Registration in progress...', (i/size(stack,3))*100);
                    Tmat(:,:,i) = eye(3,3);
                    
                    % Register current image to previous frame
                    dc = PyramidalLKOpticalFlowAffine(Mprev, squeeze(stack(:,:,i)), roi);
                    
                    % Save current image
                    Mprev = squeeze(stack(:,:,i));
                    
                    % Add current displacement to d (This is actually concatinating the two
                    % affine matrixes)
                    D = D + reshape(dc, 2, 3)*(eye(3)+[D;0 0 0]);
                    
                    % Compute displacement at current level
                    [D,k,e] = IterativeLKOpticalFlowAffine(squeeze(stack(:,:,1)), squeeze(stack(:,:,i)), roi, D);

                    % project the affine 2x3 transformation matrix onto the general affine
                    % transformation matrix
                    Tmat(1:2,:,i) = D;
                    
                    % sum up the iterations and errors
                    iter = iter + k;
                    err = err + e;
                end
            end
            
            if app.RADIO_IR_LKFlow.Value
                % this method uses a simplified Lucas-Kanade method to
                % register images solely based on their tranlatory
                % deviation. It uses a hierarchical gradient-based
                % optimization method, using 6 levels of Low-Pass filtering
                
                % create the region of continuous flow for lowest hierarchy
                % level of the gaussian pyramid
                roi=[2 2 size(stack,1)-1 size(stack,2)-1];
                
                % Tmat(:,:,1) = eye(3,3);
                
                for i=1:size(stack,3)
                    showprogress(app, 'Image Registration in progress...', (i/size(stack,3))*100);
                    Tmat(:,:,i) = eye(3,3);
                    
                    % Register current image to previous frame
                    [D, k, e] = PyramidalLKOpticalFlow(baseFrame, squeeze(stack(:,:,i)), roi);
                    
                    % project the 2D-motion vector onto the general affine
                    % transformation matrix
                    Tmat(1:2,3,i) = D.';
                    
                    % sum up the iterations and errors
                    iter = iter + k;
                    err = err + e;
                end
            end
            
            % compute the mean number of iterations (rounded) and the
            % mean relative error
            iter = ceil(iter/(size(stack,3)-1));
            err = err/(size(stack,3)-1);
        end
        
        function [HR, iter] = superResolution(app, LR, Tmat, resFactor, params, Hpsf)
            showprogress(app, ' Estimating High Resolution image', 0);
            %% Adaptive Kernel Regression
            
            if app.RADIO_SR_KernelReg.Value
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
                    showprogress(app, ' Estimating High Resolution image', (iter/maxIter*100));
                    
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
            end
            
            %% Cubic Spline Interpolation
            if app.RADIO_SR_Cubic.Value
                
                % Initialize guess as interpolated version of LR
                [X,Y]=meshgrid(0:resFactor:(size(LR,2)-1)*resFactor, 0:resFactor:(size(LR,1)-1)*resFactor);
                [XI,YI]=meshgrid(resFactor+1:(size(LR,2)-2)*resFactor-1, resFactor+1:(size(LR,1)-2)*resFactor-1);
                
                Z=interp2(X, Y, squeeze(LR(:,:,1)), XI, YI, '*spline');
                
                % Deblur the HR image and regulate using bilatural filter
                
                % Loop and improve HR in steepest descent direction
                HR = Z;
                iter = 1;
                A = ones(size(HR));
                
                %h=waitbar(0, 'Estimating high-resolution image');
                
                while iter<params.maxIter
                  showprogress(app, ' Estimating High Resolution image', (iter/params.maxIter*100));
                  %waitbar(iter/params.maxIter);
                  
                  % Compute gradient of the energy part of the cost function
                  Gback = FastGradientBackProject(HR, Z, A, Hpsf);
                
                  % Compute the gradient of the bilateral filter part of the cost function
                  Greg = GradientRegulization(HR, params.P, params.alpha);
                
                  % Perform a single SD step
                  HR = HR - params.beta.*(Gback + params.lambda.* Greg);
                  
                  iter=iter+1;
                 end
            end
            
            %% Robust SR
            if app.RADIO_SR_Robust.Value
                
                % Extract the motion vector from the affine Transformation
                % matrix
                Tvec = squeeze(Tmat(1:2, 3,:)).';
                
                % project the translation to the new image size, rounded to
                % the nearest neighbour
                Tvec_r = round(Tvec.*resFactor);
                
                %backproject the rounded vector to the initial size
                D = mod(floor(Tvec_r/resFactor),resFactor)+resFactor;
                
                % Shift all images so D is bounded from 0-resFactor
                stack = zeros(size(LR,1), size(LR,2), size(LR,3));
                [X,Y]=meshgrid(1:size(LR, 2), 1:size(LR, 1));
                
                for i=1:size(LR, 3)
                
                    stack(:,:,i)=interp2(X+Tvec_r(i,1), Y+Tvec_r(i,2), stack(:,:,i), X, Y, '*nearest');
                
                end
                
                stack_r = LR(3:end-2,3:end-2,:);
                
                % Compute initial estimate of blurred HR by the means of MedianAndShift
                [Z, A]=MedianAndShift(stack_r, D, [(size(stack_r,1)+1)*resFactor-1 (size(stack_r,2)+1)*resFactor-1], resFactor);

                % Deblur the HR image and regulate using bilatural filter
                
                % Loop and improve HR in steepest descent direction
                HR = Z;
                
                iter = 1;
                
                while iter<params.maxIter
                
                    showprogress(app, ' Estimating High Resolution image', (iter/params.maxIter*100));
%                     waitbar(iter/params.maxIter);
                    
                    % Compute gradient of the energy part of the cost function
                    Gback = GradientBackProject(HR, stack_r, D, Hpsf, resFactor);
                    
                    % Compute the gradient of the bilateral filter part of the cost function
                    Greg = GradientRegulization(HR, params.P, params.alpha);
                    
                    % Perform a single SD step
                    HR = HR - params.beta.*(Gback + params.lambda.* Greg);
                    
                    iter=iter+1;
                end
            end
            
            %% Fast Robust SR
            if app.RADIO_SR_FastRobust.Value
                
                % prepare data
                % [Z, ~, A] = APP_prepareRSR(LR, Tmat);
                
                % Extract the motion vector from the affine Transformation
                % matrix
                Tvec = squeeze(Tmat(1:2, 3,:)).';
                
                % project the translation to the new image size, rounded to
                % the nearest neighbour
                Tvec_r = round(Tvec.*resFactor);
                
                %backproject the rounded vector to the initial size
                D = mod(floor(Tvec_r/resFactor),resFactor)+resFactor;
                
                % Shift all images so D is bounded from 0-resFactor
                stack = zeros(size(LR,1), size(LR,2), size(LR,3));
                [X,Y]=meshgrid(1:size(LR, 2), 1:size(LR, 1));
                
                for i=1:size(LR, 3)
                
                    stack(:,:,i)=interp2(X+Tvec_r(i,1), Y+Tvec_r(i,2), stack(:,:,i), X, Y, '*nearest');
                
                end
                
                stack_r = LR(3:end-2,3:end-2,:);
                
                % Compute initial estimate of blurred HR by the means of MedianAndShift
                [Z, A]=MedianAndShift(stack_r, D, [(size(stack_r,1)+1)*resFactor-1 (size(stack_r,2)+1)*resFactor-1], resFactor);

                % Deblur the HR image and regulate using bilatural filter
                
                % Loop and improve HR in steepest descent direction
                HR = Z;
                iter = 1;
                
                while iter<params.maxIter
                    showprogress(app, ' Estimating High Resolution image', (iter/params.maxIter*100));
                    
                    % Compute gradient of the energy part of the cost function
                    Gback = FastGradientBackProject(HR, Z, A, Hpsf);
                    
                    % Compute the gradient of the bilateral filter part of the cost function
                    Greg = GradientRegulization(HR, params.P, params.alpha);
                    
                    % Perform a single SD step
                    HR = HR - params.beta.*(Gback + params.lambda.* Greg);
                    
                    iter=iter+1; 
                end
            end
            showprogress(app, '', 100);
        end
        
        function LR_Reg = imageRegistration(app, Tmat)
            % uses the MATLAB imwarp-function to register all images of the
            % image stack to the first frame based on a stack of affine
            % transformation matrices for each frame
            stack = app.LRstack;
            LR_Reg = zeros(size(stack,1),size(stack,2), size(stack,3));
            height = size(stack,1);
            width = size(stack,2);
            
            if(app.RADIO_IR_LKFlowAffine.Value)
                
                Tmat(1,3,2:end) = Tmat(1,3,2:end) - width/2;
                Tmat(2,3,2:end) = Tmat(2,3,2:end) + height/2;
            end
            
            for i=1:size(stack,3)
                I = 0;
                % transform Tmat into the (wrong) MATLAB transformation matrix
                trans = squeeze(Tmat(:,:,i));
                trans = trans.';
                trans(1:2,1:2) = trans(1:2,1:2).';
                
                % create the 2D transformation object needed for the
                % transformation and warp the image with an additional size restriction
                tform = affine2d(trans);
                I = imwarp(stack(:,:,i),tform,'cubic','FillValues',128);
                % LR_Reg(:,:,i) = imwarp(stack(:,:,i),tform,'cubic','OutputView', imref2d(size(stack(:,:,i)),height/2, width/2));
                ymin = floor((size(I,1)-height)/2)+1;
                xmin = floor((size(I,2)-width)/2)+1;
                LR_Reg(:,:,i) = I(ymin:ymin+height-1,xmin:xmin+width-1);
            end
            % showprogress(app, '', 100);
        end
        
        function showprogress(app, title, perc)
            app.StatusPanel.Title = strcat('Status: ', title);
            numBars = uint8(perc/10);
            
            switch(numBars)
                case 0
                    app.perc0.Visible = true;
                case 1
                    app.perc1.Visible = true;
                case 2
                    app.perc2.Visible = true;
                case 3
                    app.perc3.Visible = true;
                case 4
                    app.perc4.Visible = true;
                case 5
                    app.perc5.Visible = true;
                case 6
                    app.perc6.Visible = true;
                case 7
                    app.perc7.Visible = true;
                case 8
                    app.perc8.Visible = true;
                case 9
                    app.perc9.Visible = true;
                case 10
                    app.perc0.Visible = false;
                    app.perc1.Visible = false;
                    app.perc2.Visible = false;
                    app.perc3.Visible = false;
                    app.perc4.Visible = false;
                    app.perc5.Visible = false;
                    app.perc6.Visible = false;
                    app.perc7.Visible = false;
                    app.perc8.Visible = false;
                    app.perc9.Visible = false;
                    app.StatusPanel.Title = 'Status: ';
            end
            drawnow
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % Load all components
            addpath([pwd '/Helper']);
            % Image Registration
            addpath([pwd '/ImageRegistration/LKOFlow']);
            addpath([pwd '/ImageRegistration/LKOFlowAffine']);
            % Super Resolution
            addpath([pwd '/SuperResolution/SplineInterpolation']);
            addpath([pwd '/SuperResolution/AdaptiveKernel']);
            addpath([pwd '/SuperResolution/Robust']);
            addpath([pwd '/SuperResolution/FastRobust']);
            
            app.imReg_flag = false;
        end

        % Button pushed function: BTN_loadFile
        function BTN_loadFileButtonPushed(app, event)
            filter = '*.avi;*.mov;*.mp4;.m4v;';
            [FileName,PathName] = uigetfile(filter,'Select the movie file (avi, mov, mp4, m4v)');
            drawnow;
            figure(app.MFSRToolUIFigure);
            
            if FileName ~= 0
                % Load video file
                inVideo = LoadVideo([PathName FileName]);
                app.LRstack = inVideo;
                app.frameCnt = 1;
                
                app.IMG_LR.ImageSource = APP_buildIMG(inVideo(:,:,1));
                
                % Fill the sequence info table
                app.VAL_imInfo_format.Text = FileName(end-3:end);
                app.VAL_imInfo_nframes.Text = int2str(size(inVideo,3));
                app.VAL_imInfo_imHeight.Text = int2str(size(inVideo,1));
                app.VAL_imInfo_imWidth.Text = int2str(size(inVideo,2));
                
                % Enable buttons
                app.BTN_next.Enable = 'on';
                app.BTN_enhance.Enable = 'on';
                app.BTN_saveLRFrame.Enable = 'on';
                
                % Reset the IR-Flag
                app.imReg_flag = false;
            end
        end

        % Button pushed function: BTN_next
        function BTN_nextButtonPushed(app, event)
            
            % Display the next frame in the stack
            app.frameCnt = app.frameCnt + 1;
            app.IMG_LR.ImageSource = APP_buildIMG(app.LRstack(:,:,app.frameCnt));
            
            % If we hit the ceiling, disable 'next' button
            if app.frameCnt == size(app.LRstack,3)
                app.BTN_next.Enable = 'off';
            end
            
            % if we're past the first frame, enable 'previus' button
            if app.frameCnt > 1
                app.BTN_prev.Enable = 'on';
            end
        end

        % Button pushed function: BTN_prev
        function BTN_prevButtonPushed(app, event)
            
            % Display the previous frame in the stack
            app.frameCnt = app.frameCnt - 1;
            app.IMG_LR.ImageSource = APP_buildIMG(app.LRstack(:,:,app.frameCnt));
            
            % if we hit the bottom, disable 'previous' button
            if app.frameCnt == 1
                app.BTN_prev.Enable = 'off';
            end
            
            % if we're below the last frame, enable 'next' button
            if app.frameCnt < size(app.LRstack,3)
                app.BTN_next.Enable = 'on';
            end
        end

        % Button pushed function: BTN_saveLRFrame
        function BTN_saveLRFrameButtonPushed(app, event)
            [FileName,PathName] = uiputfile('*.jpg','Save image file');
            drawnow;
            figure(app.MFSRToolUIFigure);
            
            if FileName ~= 0
              imwrite(uint8(app.LRstack(:,:,app.frameCnt)), [PathName FileName]);
            end
        end

        % Button pushed function: BTN_enhance
        function BTN_enhanceButtonPushed(app, event)
            
            % Collect the parameters from the GUI
            [params, resFactor, Hpsf] = collectParams(app);
            
            % Compute the affine transformation matrices for each frame of
            % the sequence using the method chosen in the GUI
            if(~app.imReg_flag)
                tic;
                [app.TM, iter, err] = getTransform(app);
                
                % fill the benchmark measurement table
                app.VAL_IR_t.Text = num2str(toc);
                app.VAL_IR_err.Text = num2str(err);
                app.VAL_IR_n.Text = num2str(iter);
                
                % Actually execute the image registration using the stack of
                % geometric transformations in TM
                app.LR_reg = imageRegistration(app,app.TM);
                for i=1:size(app.LR_reg,3)
                    figure(i)
                    debug = app.LR_reg(:,:,i);
                    image(debug)
                end

                app.imReg_flag = true;
            end
            
            tic;
            [app.HRimage, iter] = superResolution(app, app.LR_reg, app.TM, resFactor, params, Hpsf);
            
            % fill the benchmark measurement table
            app.VAL_SR_t.Text = num2str(toc);
            app.VAL_SR_err.Text = 'n.A.'; % num2str(err);
            app.VAL_SR_n.Text = num2str(iter);
            
            app.IMG_HR.ImageSource = APP_buildIMG(app.HRimage);
            
            app.BTN_saveHRFrame.Enable = true;
            
        end

        % Button pushed function: BTN_saveHRFrame
        function BTN_saveHRFrameButtonPushed(app, event)
            [FileName,PathName] = uiputfile('*.jpg','Save image file');
            drawnow;
            figure(app.MFSRToolUIFigure);
            
            if FileName ~= 0
              imwrite(uint8(app.HRimage), [PathName FileName]);
            end
        end

        % Selection changed function: RADIOGROUP_IR_Method
        function RADIOGROUP_IR_MethodSelectionChanged(app, event)
            app.imReg_flag = false;
            
        end

        % Button pushed function: BTN_reset
        function BTN_resetButtonPushed(app, event)
            app.imReg_flag = false;
            app.IMG_HR.ImageSource = '';
            drawnow
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create MFSRToolUIFigure and hide until all components are created
            app.MFSRToolUIFigure = uifigure('Visible', 'off');
            app.MFSRToolUIFigure.Color = [0.149 0.149 0.149];
            app.MFSRToolUIFigure.Position = [100 100 800 800];
            app.MFSRToolUIFigure.Name = 'MFSR Tool';

            % Create IMG_LR
            app.IMG_LR = uiimage(app.MFSRToolUIFigure);
            app.IMG_LR.Position = [30 400 350 350];

            % Create IMG_HR
            app.IMG_HR = uiimage(app.MFSRToolUIFigure);
            app.IMG_HR.Position = [420 400 350 350];

            % Create LABEL_heading_LRstack
            app.LABEL_heading_LRstack = uilabel(app.MFSRToolUIFigure);
            app.LABEL_heading_LRstack.FontSize = 20;
            app.LABEL_heading_LRstack.FontWeight = 'bold';
            app.LABEL_heading_LRstack.FontColor = [1 1 0];
            app.LABEL_heading_LRstack.Position = [31 757 214 32];
            app.LABEL_heading_LRstack.Text = 'Low Resolution Stack';

            % Create LABEL_heading_HRstack
            app.LABEL_heading_HRstack = uilabel(app.MFSRToolUIFigure);
            app.LABEL_heading_HRstack.FontSize = 20;
            app.LABEL_heading_HRstack.FontWeight = 'bold';
            app.LABEL_heading_HRstack.FontColor = [1 1 0];
            app.LABEL_heading_HRstack.Position = [431 757 240 32];
            app.LABEL_heading_HRstack.Text = 'Super Resolution Result';

            % Create PANEL_imInfo
            app.PANEL_imInfo = uipanel(app.MFSRToolUIFigure);
            app.PANEL_imInfo.ForegroundColor = [1 1 0];
            app.PANEL_imInfo.Title = 'Low Res Image Information';
            app.PANEL_imInfo.BackgroundColor = [0 0 0];
            app.PANEL_imInfo.FontWeight = 'bold';
            app.PANEL_imInfo.FontSize = 16;
            app.PANEL_imInfo.Position = [30 230 350 120];

            % Create GRID_imInfo
            app.GRID_imInfo = uigridlayout(app.PANEL_imInfo);
            app.GRID_imInfo.ColumnWidth = {'1x', '1x', '1x', '1x'};

            % Create LABEL_imInfo_nframes
            app.LABEL_imInfo_nframes = uilabel(app.GRID_imInfo);
            app.LABEL_imInfo_nframes.BackgroundColor = [0 0 0];
            app.LABEL_imInfo_nframes.FontSize = 14;
            app.LABEL_imInfo_nframes.FontWeight = 'bold';
            app.LABEL_imInfo_nframes.FontColor = [1 1 1];
            app.LABEL_imInfo_nframes.Layout.Row = 1;
            app.LABEL_imInfo_nframes.Layout.Column = 1;
            app.LABEL_imInfo_nframes.Text = 'Frames:';

            % Create LABEL_imInfo_imHeight
            app.LABEL_imInfo_imHeight = uilabel(app.GRID_imInfo);
            app.LABEL_imInfo_imHeight.BackgroundColor = [0 0 0];
            app.LABEL_imInfo_imHeight.FontSize = 14;
            app.LABEL_imInfo_imHeight.FontWeight = 'bold';
            app.LABEL_imInfo_imHeight.FontColor = [1 1 1];
            app.LABEL_imInfo_imHeight.Layout.Row = 2;
            app.LABEL_imInfo_imHeight.Layout.Column = 3;
            app.LABEL_imInfo_imHeight.Text = 'Height (px):';

            % Create LABEL_imInfo_format
            app.LABEL_imInfo_format = uilabel(app.GRID_imInfo);
            app.LABEL_imInfo_format.BackgroundColor = [0 0 0];
            app.LABEL_imInfo_format.FontSize = 14;
            app.LABEL_imInfo_format.FontWeight = 'bold';
            app.LABEL_imInfo_format.FontColor = [1 1 1];
            app.LABEL_imInfo_format.Layout.Row = 1;
            app.LABEL_imInfo_format.Layout.Column = 3;
            app.LABEL_imInfo_format.Text = 'Format:';

            % Create LABEL_imInfo_imWidth
            app.LABEL_imInfo_imWidth = uilabel(app.GRID_imInfo);
            app.LABEL_imInfo_imWidth.BackgroundColor = [0 0 0];
            app.LABEL_imInfo_imWidth.FontSize = 14;
            app.LABEL_imInfo_imWidth.FontWeight = 'bold';
            app.LABEL_imInfo_imWidth.FontColor = [1 1 1];
            app.LABEL_imInfo_imWidth.Layout.Row = 2;
            app.LABEL_imInfo_imWidth.Layout.Column = 1;
            app.LABEL_imInfo_imWidth.Text = 'Width (px):';

            % Create VAL_imInfo_nframes
            app.VAL_imInfo_nframes = uilabel(app.GRID_imInfo);
            app.VAL_imInfo_nframes.BackgroundColor = [0.149 0.149 0.149];
            app.VAL_imInfo_nframes.HorizontalAlignment = 'right';
            app.VAL_imInfo_nframes.FontSize = 14;
            app.VAL_imInfo_nframes.FontColor = [0.302 0.7451 0.9333];
            app.VAL_imInfo_nframes.Layout.Row = 1;
            app.VAL_imInfo_nframes.Layout.Column = 2;
            app.VAL_imInfo_nframes.Text = '';

            % Create VAL_imInfo_format
            app.VAL_imInfo_format = uilabel(app.GRID_imInfo);
            app.VAL_imInfo_format.BackgroundColor = [0.149 0.149 0.149];
            app.VAL_imInfo_format.HorizontalAlignment = 'right';
            app.VAL_imInfo_format.FontSize = 14;
            app.VAL_imInfo_format.FontColor = [0.302 0.7451 0.9333];
            app.VAL_imInfo_format.Layout.Row = 1;
            app.VAL_imInfo_format.Layout.Column = 4;
            app.VAL_imInfo_format.Text = '';

            % Create VAL_imInfo_imWidth
            app.VAL_imInfo_imWidth = uilabel(app.GRID_imInfo);
            app.VAL_imInfo_imWidth.BackgroundColor = [0.149 0.149 0.149];
            app.VAL_imInfo_imWidth.HorizontalAlignment = 'right';
            app.VAL_imInfo_imWidth.FontSize = 14;
            app.VAL_imInfo_imWidth.FontColor = [0.302 0.7451 0.9333];
            app.VAL_imInfo_imWidth.Layout.Row = 2;
            app.VAL_imInfo_imWidth.Layout.Column = 2;
            app.VAL_imInfo_imWidth.Text = '';

            % Create VAL_imInfo_imHeight
            app.VAL_imInfo_imHeight = uilabel(app.GRID_imInfo);
            app.VAL_imInfo_imHeight.BackgroundColor = [0.149 0.149 0.149];
            app.VAL_imInfo_imHeight.HorizontalAlignment = 'right';
            app.VAL_imInfo_imHeight.FontSize = 14;
            app.VAL_imInfo_imHeight.FontColor = [0.302 0.7451 0.9333];
            app.VAL_imInfo_imHeight.Layout.Row = 2;
            app.VAL_imInfo_imHeight.Layout.Column = 4;
            app.VAL_imInfo_imHeight.Text = '';

            % Create BTN_next
            app.BTN_next = uibutton(app.MFSRToolUIFigure, 'push');
            app.BTN_next.ButtonPushedFcn = createCallbackFcn(app, @BTN_nextButtonPushed, true);
            app.BTN_next.BackgroundColor = [0 0 0];
            app.BTN_next.FontSize = 14;
            app.BTN_next.FontWeight = 'bold';
            app.BTN_next.FontColor = [1 1 0];
            app.BTN_next.Enable = 'off';
            app.BTN_next.Position = [321 361 50 30];
            app.BTN_next.Text = '>>';

            % Create BTN_prev
            app.BTN_prev = uibutton(app.MFSRToolUIFigure, 'push');
            app.BTN_prev.ButtonPushedFcn = createCallbackFcn(app, @BTN_prevButtonPushed, true);
            app.BTN_prev.BackgroundColor = [0 0 0];
            app.BTN_prev.FontSize = 14;
            app.BTN_prev.FontWeight = 'bold';
            app.BTN_prev.FontColor = [1 1 0];
            app.BTN_prev.Enable = 'off';
            app.BTN_prev.Position = [31 361 50 30];
            app.BTN_prev.Text = '<<';

            % Create BTN_loadFile
            app.BTN_loadFile = uibutton(app.MFSRToolUIFigure, 'push');
            app.BTN_loadFile.ButtonPushedFcn = createCallbackFcn(app, @BTN_loadFileButtonPushed, true);
            app.BTN_loadFile.BackgroundColor = [0 0 0];
            app.BTN_loadFile.FontSize = 14;
            app.BTN_loadFile.FontWeight = 'bold';
            app.BTN_loadFile.FontColor = [1 1 1];
            app.BTN_loadFile.Position = [91 361 120 30];
            app.BTN_loadFile.Text = 'Load Video';

            % Create PANEL_bench
            app.PANEL_bench = uipanel(app.MFSRToolUIFigure);
            app.PANEL_bench.ForegroundColor = [1 1 0];
            app.PANEL_bench.Title = 'Benchmark';
            app.PANEL_bench.BackgroundColor = [0 0 0];
            app.PANEL_bench.FontWeight = 'bold';
            app.PANEL_bench.FontSize = 16;
            app.PANEL_bench.Position = [31 21 350 200];

            % Create GRID_bench
            app.GRID_bench = uigridlayout(app.PANEL_bench);
            app.GRID_bench.RowHeight = {'1x', '4x'};

            % Create LABEL_heading_IR
            app.LABEL_heading_IR = uilabel(app.GRID_bench);
            app.LABEL_heading_IR.BackgroundColor = [0 0 0];
            app.LABEL_heading_IR.VerticalAlignment = 'bottom';
            app.LABEL_heading_IR.FontSize = 14;
            app.LABEL_heading_IR.FontWeight = 'bold';
            app.LABEL_heading_IR.FontColor = [0.302 0.7451 0.9333];
            app.LABEL_heading_IR.Layout.Row = 1;
            app.LABEL_heading_IR.Layout.Column = 1;
            app.LABEL_heading_IR.Text = 'Image Registration';

            % Create LABEL_heading_SR
            app.LABEL_heading_SR = uilabel(app.GRID_bench);
            app.LABEL_heading_SR.BackgroundColor = [0 0 0];
            app.LABEL_heading_SR.VerticalAlignment = 'bottom';
            app.LABEL_heading_SR.FontSize = 14;
            app.LABEL_heading_SR.FontWeight = 'bold';
            app.LABEL_heading_SR.FontColor = [0.302 0.7451 0.9333];
            app.LABEL_heading_SR.Layout.Row = 1;
            app.LABEL_heading_SR.Layout.Column = 2;
            app.LABEL_heading_SR.Text = 'Super Resolution';

            % Create GRID_IR_bench
            app.GRID_IR_bench = uigridlayout(app.GRID_bench);
            app.GRID_IR_bench.RowHeight = {'1x', '1x', '1x'};
            app.GRID_IR_bench.Layout.Row = 2;
            app.GRID_IR_bench.Layout.Column = 1;

            % Create LABEL_IR_time
            app.LABEL_IR_time = uilabel(app.GRID_IR_bench);
            app.LABEL_IR_time.BackgroundColor = [0 0 0];
            app.LABEL_IR_time.FontWeight = 'bold';
            app.LABEL_IR_time.FontColor = [1 1 1];
            app.LABEL_IR_time.Layout.Row = 1;
            app.LABEL_IR_time.Layout.Column = 1;
            app.LABEL_IR_time.Text = 'time (sec):';

            % Create LABEL_IR_iter
            app.LABEL_IR_iter = uilabel(app.GRID_IR_bench);
            app.LABEL_IR_iter.BackgroundColor = [0 0 0];
            app.LABEL_IR_iter.FontWeight = 'bold';
            app.LABEL_IR_iter.FontColor = [1 1 1];
            app.LABEL_IR_iter.Layout.Row = 2;
            app.LABEL_IR_iter.Layout.Column = 1;
            app.LABEL_IR_iter.Text = 'iterations:';

            % Create LABEL_IR_err
            app.LABEL_IR_err = uilabel(app.GRID_IR_bench);
            app.LABEL_IR_err.BackgroundColor = [0 0 0];
            app.LABEL_IR_err.FontWeight = 'bold';
            app.LABEL_IR_err.FontColor = [1 1 1];
            app.LABEL_IR_err.Layout.Row = 3;
            app.LABEL_IR_err.Layout.Column = 1;
            app.LABEL_IR_err.Text = 'error:';

            % Create VAL_IR_t
            app.VAL_IR_t = uilabel(app.GRID_IR_bench);
            app.VAL_IR_t.BackgroundColor = [0.149 0.149 0.149];
            app.VAL_IR_t.HorizontalAlignment = 'right';
            app.VAL_IR_t.FontAngle = 'italic';
            app.VAL_IR_t.FontColor = [0.302 0.7451 0.9333];
            app.VAL_IR_t.Layout.Row = 1;
            app.VAL_IR_t.Layout.Column = 2;
            app.VAL_IR_t.Text = '';

            % Create VAL_IR_n
            app.VAL_IR_n = uilabel(app.GRID_IR_bench);
            app.VAL_IR_n.BackgroundColor = [0.149 0.149 0.149];
            app.VAL_IR_n.HorizontalAlignment = 'right';
            app.VAL_IR_n.FontAngle = 'italic';
            app.VAL_IR_n.FontColor = [0.302 0.7451 0.9333];
            app.VAL_IR_n.Layout.Row = 2;
            app.VAL_IR_n.Layout.Column = 2;
            app.VAL_IR_n.Text = '';

            % Create VAL_IR_err
            app.VAL_IR_err = uilabel(app.GRID_IR_bench);
            app.VAL_IR_err.BackgroundColor = [0.149 0.149 0.149];
            app.VAL_IR_err.HorizontalAlignment = 'right';
            app.VAL_IR_err.FontAngle = 'italic';
            app.VAL_IR_err.FontColor = [0.302 0.7451 0.9333];
            app.VAL_IR_err.Layout.Row = 3;
            app.VAL_IR_err.Layout.Column = 2;
            app.VAL_IR_err.Text = '';

            % Create GRID_SR_bench
            app.GRID_SR_bench = uigridlayout(app.GRID_bench);
            app.GRID_SR_bench.RowHeight = {'1x', '1x', '1x'};
            app.GRID_SR_bench.Layout.Row = 2;
            app.GRID_SR_bench.Layout.Column = 2;

            % Create LABEL_SR_time
            app.LABEL_SR_time = uilabel(app.GRID_SR_bench);
            app.LABEL_SR_time.BackgroundColor = [0 0 0];
            app.LABEL_SR_time.FontWeight = 'bold';
            app.LABEL_SR_time.FontColor = [1 1 1];
            app.LABEL_SR_time.Layout.Row = 1;
            app.LABEL_SR_time.Layout.Column = 1;
            app.LABEL_SR_time.Text = 'time (ms):';

            % Create LABEL_SR_iter
            app.LABEL_SR_iter = uilabel(app.GRID_SR_bench);
            app.LABEL_SR_iter.BackgroundColor = [0 0 0];
            app.LABEL_SR_iter.FontWeight = 'bold';
            app.LABEL_SR_iter.FontColor = [1 1 1];
            app.LABEL_SR_iter.Layout.Row = 2;
            app.LABEL_SR_iter.Layout.Column = 1;
            app.LABEL_SR_iter.Text = 'iterations:';

            % Create LABEL_SR_err
            app.LABEL_SR_err = uilabel(app.GRID_SR_bench);
            app.LABEL_SR_err.BackgroundColor = [0 0 0];
            app.LABEL_SR_err.FontWeight = 'bold';
            app.LABEL_SR_err.FontColor = [1 1 1];
            app.LABEL_SR_err.Layout.Row = 3;
            app.LABEL_SR_err.Layout.Column = 1;
            app.LABEL_SR_err.Text = 'error';

            % Create VAL_SR_t
            app.VAL_SR_t = uilabel(app.GRID_SR_bench);
            app.VAL_SR_t.BackgroundColor = [0.149 0.149 0.149];
            app.VAL_SR_t.HorizontalAlignment = 'right';
            app.VAL_SR_t.FontAngle = 'italic';
            app.VAL_SR_t.FontColor = [0.302 0.7451 0.9333];
            app.VAL_SR_t.Layout.Row = 1;
            app.VAL_SR_t.Layout.Column = 2;
            app.VAL_SR_t.Text = '';

            % Create VAL_SR_n
            app.VAL_SR_n = uilabel(app.GRID_SR_bench);
            app.VAL_SR_n.BackgroundColor = [0.149 0.149 0.149];
            app.VAL_SR_n.HorizontalAlignment = 'right';
            app.VAL_SR_n.FontAngle = 'italic';
            app.VAL_SR_n.FontColor = [0.302 0.7451 0.9333];
            app.VAL_SR_n.Layout.Row = 2;
            app.VAL_SR_n.Layout.Column = 2;
            app.VAL_SR_n.Text = '';

            % Create VAL_SR_err
            app.VAL_SR_err = uilabel(app.GRID_SR_bench);
            app.VAL_SR_err.BackgroundColor = [0.149 0.149 0.149];
            app.VAL_SR_err.HorizontalAlignment = 'right';
            app.VAL_SR_err.FontAngle = 'italic';
            app.VAL_SR_err.FontColor = [0.302 0.7451 0.9333];
            app.VAL_SR_err.Layout.Row = 3;
            app.VAL_SR_err.Layout.Column = 2;
            app.VAL_SR_err.Text = '';

            % Create BTN_enhance
            app.BTN_enhance = uibutton(app.MFSRToolUIFigure, 'push');
            app.BTN_enhance.ButtonPushedFcn = createCallbackFcn(app, @BTN_enhanceButtonPushed, true);
            app.BTN_enhance.BackgroundColor = [0 0 0];
            app.BTN_enhance.FontSize = 14;
            app.BTN_enhance.FontWeight = 'bold';
            app.BTN_enhance.FontColor = [1 1 1];
            app.BTN_enhance.Enable = 'off';
            app.BTN_enhance.Position = [521 361 130 30];
            app.BTN_enhance.Text = 'ENHANCE';

            % Create RADIOGROUP_IR_Method
            app.RADIOGROUP_IR_Method = uibuttongroup(app.MFSRToolUIFigure);
            app.RADIOGROUP_IR_Method.SelectionChangedFcn = createCallbackFcn(app, @RADIOGROUP_IR_MethodSelectionChanged, true);
            app.RADIOGROUP_IR_Method.ForegroundColor = [1 1 0];
            app.RADIOGROUP_IR_Method.Title = 'Registration Method';
            app.RADIOGROUP_IR_Method.BackgroundColor = [0 0 0];
            app.RADIOGROUP_IR_Method.FontWeight = 'bold';
            app.RADIOGROUP_IR_Method.FontSize = 16;
            app.RADIOGROUP_IR_Method.Position = [421 230 180 121];

            % Create RADIO_IR_MATLAB
            app.RADIO_IR_MATLAB = uiradiobutton(app.RADIOGROUP_IR_Method);
            app.RADIO_IR_MATLAB.Text = 'Standard MATLAB';
            app.RADIO_IR_MATLAB.FontSize = 14;
            app.RADIO_IR_MATLAB.FontColor = [0.9412 0.9412 0.9412];
            app.RADIO_IR_MATLAB.Position = [11 70 138 22];
            app.RADIO_IR_MATLAB.Value = true;

            % Create RADIO_IR_LKFlowAffine
            app.RADIO_IR_LKFlowAffine = uiradiobutton(app.RADIOGROUP_IR_Method);
            app.RADIO_IR_LKFlowAffine.Text = 'Lucas-Kanade Affine';
            app.RADIO_IR_LKFlowAffine.FontSize = 14;
            app.RADIO_IR_LKFlowAffine.FontColor = [0.9412 0.9412 0.9412];
            app.RADIO_IR_LKFlowAffine.Position = [11 43 155 22];

            % Create RADIO_IR_LKFlow
            app.RADIO_IR_LKFlow = uiradiobutton(app.RADIOGROUP_IR_Method);
            app.RADIO_IR_LKFlow.Text = 'Lucas-Kanade Motion';
            app.RADIO_IR_LKFlow.FontSize = 14;
            app.RADIO_IR_LKFlow.FontColor = [0.9412 0.9412 0.9412];
            app.RADIO_IR_LKFlow.Position = [11 16 158 22];

            % Create RADIOGROUP_SR_Method
            app.RADIOGROUP_SR_Method = uibuttongroup(app.MFSRToolUIFigure);
            app.RADIOGROUP_SR_Method.ForegroundColor = [1 1 0];
            app.RADIOGROUP_SR_Method.Title = 'Super Resolution Method';
            app.RADIOGROUP_SR_Method.BackgroundColor = [0 0 0];
            app.RADIOGROUP_SR_Method.FontWeight = 'bold';
            app.RADIOGROUP_SR_Method.FontSize = 16;
            app.RADIOGROUP_SR_Method.Position = [421 81 350 140];

            % Create RADIO_SR_KernelReg
            app.RADIO_SR_KernelReg = uiradiobutton(app.RADIOGROUP_SR_Method);
            app.RADIO_SR_KernelReg.Text = 'Adaptive Kernel Regression';
            app.RADIO_SR_KernelReg.FontSize = 14;
            app.RADIO_SR_KernelReg.FontColor = [0.9412 0.9412 0.9412];
            app.RADIO_SR_KernelReg.Position = [11 89 197 22];
            app.RADIO_SR_KernelReg.Value = true;

            % Create RADIO_SR_Cubic
            app.RADIO_SR_Cubic = uiradiobutton(app.RADIOGROUP_SR_Method);
            app.RADIO_SR_Cubic.Text = 'Cubic Spline Interpolation';
            app.RADIO_SR_Cubic.FontSize = 14;
            app.RADIO_SR_Cubic.FontColor = [0.9412 0.9412 0.9412];
            app.RADIO_SR_Cubic.Position = [11 62 182 22];

            % Create RADIO_SR_Robust
            app.RADIO_SR_Robust = uiradiobutton(app.RADIOGROUP_SR_Method);
            app.RADIO_SR_Robust.Text = 'Robust Super Resolution';
            app.RADIO_SR_Robust.FontSize = 14;
            app.RADIO_SR_Robust.FontColor = [0.9412 0.9412 0.9412];
            app.RADIO_SR_Robust.Position = [11 35 179 22];

            % Create RADIO_SR_FastRobust
            app.RADIO_SR_FastRobust = uiradiobutton(app.RADIOGROUP_SR_Method);
            app.RADIO_SR_FastRobust.Text = 'Fast Robust Super Resolution';
            app.RADIO_SR_FastRobust.FontSize = 14;
            app.RADIO_SR_FastRobust.FontColor = [0.9412 0.9412 0.9412];
            app.RADIO_SR_FastRobust.Position = [11 8 210 22];

            % Create BTN_reset
            app.BTN_reset = uibutton(app.MFSRToolUIFigure, 'push');
            app.BTN_reset.ButtonPushedFcn = createCallbackFcn(app, @BTN_resetButtonPushed, true);
            app.BTN_reset.BackgroundColor = [0 0 0];
            app.BTN_reset.FontSize = 14;
            app.BTN_reset.FontWeight = 'bold';
            app.BTN_reset.FontColor = [0.9294 0.6941 0.1255];
            app.BTN_reset.Position = [421 361 70 30];
            app.BTN_reset.Text = 'Reset';

            % Create PANEL_SR_params
            app.PANEL_SR_params = uipanel(app.MFSRToolUIFigure);
            app.PANEL_SR_params.ForegroundColor = [1 1 0];
            app.PANEL_SR_params.Title = 'Parameters';
            app.PANEL_SR_params.BackgroundColor = [0 0 0];
            app.PANEL_SR_params.FontWeight = 'bold';
            app.PANEL_SR_params.FontSize = 16;
            app.PANEL_SR_params.Position = [611 232 160 119];

            % Create ResfactorLabel
            app.ResfactorLabel = uilabel(app.PANEL_SR_params);
            app.ResfactorLabel.FontSize = 14;
            app.ResfactorLabel.FontColor = [0.9412 0.9412 0.9412];
            app.ResfactorLabel.Position = [10 57 90 22];
            app.ResfactorLabel.Text = 'Res - factor';

            % Create PARAM_ResFactor
            app.PARAM_ResFactor = uieditfield(app.PANEL_SR_params, 'numeric');
            app.PARAM_ResFactor.Limits = [1 10];
            app.PARAM_ResFactor.RoundFractionalValues = 'on';
            app.PARAM_ResFactor.ValueDisplayFormat = '%.0f';
            app.PARAM_ResFactor.FontSize = 14;
            app.PARAM_ResFactor.FontColor = [0.302 0.7451 0.9333];
            app.PARAM_ResFactor.BackgroundColor = [0.149 0.149 0.149];
            app.PARAM_ResFactor.Position = [120 57 31 22];
            app.PARAM_ResFactor.Value = 2;

            % Create IterationsEditFieldLabel
            app.IterationsEditFieldLabel = uilabel(app.PANEL_SR_params);
            app.IterationsEditFieldLabel.FontSize = 14;
            app.IterationsEditFieldLabel.FontColor = [0.9412 0.9412 0.9412];
            app.IterationsEditFieldLabel.Position = [11 22 90 22];
            app.IterationsEditFieldLabel.Text = 'Iterations';

            % Create PARAM_Iterations
            app.PARAM_Iterations = uieditfield(app.PANEL_SR_params, 'numeric');
            app.PARAM_Iterations.Limits = [1 9999];
            app.PARAM_Iterations.RoundFractionalValues = 'on';
            app.PARAM_Iterations.ValueDisplayFormat = '%.0f';
            app.PARAM_Iterations.FontSize = 14;
            app.PARAM_Iterations.FontColor = [0.302 0.7451 0.9333];
            app.PARAM_Iterations.BackgroundColor = [0.149 0.149 0.149];
            app.PARAM_Iterations.Position = [100 22 51 22];
            app.PARAM_Iterations.Value = 20;

            % Create BTN_saveLRFrame
            app.BTN_saveLRFrame = uibutton(app.MFSRToolUIFigure, 'push');
            app.BTN_saveLRFrame.ButtonPushedFcn = createCallbackFcn(app, @BTN_saveLRFrameButtonPushed, true);
            app.BTN_saveLRFrame.BackgroundColor = [0 0 0];
            app.BTN_saveLRFrame.FontSize = 14;
            app.BTN_saveLRFrame.FontWeight = 'bold';
            app.BTN_saveLRFrame.FontColor = [0.302 0.7451 0.9333];
            app.BTN_saveLRFrame.Enable = 'off';
            app.BTN_saveLRFrame.Position = [221 361 90 30];
            app.BTN_saveLRFrame.Text = 'Save Frame';

            % Create BTN_saveHRFrame
            app.BTN_saveHRFrame = uibutton(app.MFSRToolUIFigure, 'push');
            app.BTN_saveHRFrame.ButtonPushedFcn = createCallbackFcn(app, @BTN_saveHRFrameButtonPushed, true);
            app.BTN_saveHRFrame.BackgroundColor = [0 0 0];
            app.BTN_saveHRFrame.FontSize = 14;
            app.BTN_saveHRFrame.FontWeight = 'bold';
            app.BTN_saveHRFrame.FontColor = [0.302 0.7451 0.9333];
            app.BTN_saveHRFrame.Enable = 'off';
            app.BTN_saveHRFrame.Position = [681 361 94 30];
            app.BTN_saveHRFrame.Text = 'Save Frame';

            % Create StatusPanel
            app.StatusPanel = uipanel(app.MFSRToolUIFigure);
            app.StatusPanel.ForegroundColor = [1 1 0];
            app.StatusPanel.BorderType = 'none';
            app.StatusPanel.Title = 'Status';
            app.StatusPanel.BackgroundColor = [0.149 0.149 0.149];
            app.StatusPanel.Position = [421 21 350 50];

            % Create GridLayout
            app.GridLayout = uigridlayout(app.StatusPanel);
            app.GridLayout.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout.RowHeight = {'1x'};

            % Create perc0
            app.perc0 = uilabel(app.GridLayout);
            app.perc0.BackgroundColor = [0.302 0.7451 0.9333];
            app.perc0.Visible = 'off';
            app.perc0.Layout.Row = 1;
            app.perc0.Layout.Column = 1;
            app.perc0.Text = '';

            % Create perc1
            app.perc1 = uilabel(app.GridLayout);
            app.perc1.BackgroundColor = [0.302 0.7451 0.9333];
            app.perc1.Visible = 'off';
            app.perc1.Layout.Row = 1;
            app.perc1.Layout.Column = 2;
            app.perc1.Text = '';

            % Create perc2
            app.perc2 = uilabel(app.GridLayout);
            app.perc2.BackgroundColor = [0.302 0.7451 0.9333];
            app.perc2.Visible = 'off';
            app.perc2.Layout.Row = 1;
            app.perc2.Layout.Column = 3;
            app.perc2.Text = '';

            % Create perc3
            app.perc3 = uilabel(app.GridLayout);
            app.perc3.BackgroundColor = [0.302 0.7451 0.9333];
            app.perc3.Visible = 'off';
            app.perc3.Layout.Row = 1;
            app.perc3.Layout.Column = 4;
            app.perc3.Text = '';

            % Create perc4
            app.perc4 = uilabel(app.GridLayout);
            app.perc4.BackgroundColor = [0.302 0.7451 0.9333];
            app.perc4.Visible = 'off';
            app.perc4.Layout.Row = 1;
            app.perc4.Layout.Column = 5;
            app.perc4.Text = '';

            % Create perc5
            app.perc5 = uilabel(app.GridLayout);
            app.perc5.BackgroundColor = [0.302 0.7451 0.9333];
            app.perc5.Visible = 'off';
            app.perc5.Layout.Row = 1;
            app.perc5.Layout.Column = 6;
            app.perc5.Text = '';

            % Create perc6
            app.perc6 = uilabel(app.GridLayout);
            app.perc6.BackgroundColor = [0.302 0.7451 0.9333];
            app.perc6.Visible = 'off';
            app.perc6.Layout.Row = 1;
            app.perc6.Layout.Column = 7;
            app.perc6.Text = '';

            % Create perc7
            app.perc7 = uilabel(app.GridLayout);
            app.perc7.BackgroundColor = [0.302 0.7451 0.9333];
            app.perc7.Visible = 'off';
            app.perc7.Layout.Row = 1;
            app.perc7.Layout.Column = 8;
            app.perc7.Text = '';

            % Create perc8
            app.perc8 = uilabel(app.GridLayout);
            app.perc8.BackgroundColor = [0.302 0.7451 0.9333];
            app.perc8.Visible = 'off';
            app.perc8.Layout.Row = 1;
            app.perc8.Layout.Column = 9;
            app.perc8.Text = '';

            % Create perc9
            app.perc9 = uilabel(app.GridLayout);
            app.perc9.BackgroundColor = [0.302 0.7451 0.9333];
            app.perc9.Visible = 'off';
            app.perc9.Layout.Row = 1;
            app.perc9.Layout.Column = 10;
            app.perc9.Text = '';

            % Show the figure after all components are created
            app.MFSRToolUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = MFSR_App

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.MFSRToolUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.MFSRToolUIFigure)
        end
    end
end