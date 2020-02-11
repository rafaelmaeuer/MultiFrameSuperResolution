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
        ResFactorLabel            matlab.ui.control.Label
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
        
        % Get parameters from GUI and set some internal
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
        
        % Returns the translation vector for each frame f with respect to the
        % reference frame (first frame in stack) and the registered stack of images
        function [LR_reg, Tvec, iter, err] = getTransform(app)
            
            % initialize Progress bar
            ShowProgress(app, 'Image Registration in progress...', 0);
            
            % load the image stack handle to a local variable to increase speed
            stack = app.LRstack;
            
            % initialize output parameters
            LR_reg = zeros(size(stack));
            Tvec = zeros(size(stack,3),2);
            iter = 0;
            err = 0;
            
            %% Matlab Image Registration
            if app.RADIO_IR_MATLAB.Value
                [LR_reg, Tvec] = RegisterImageSeqMatlab(app, stack);
            end
            
            %% LK Optical Flow Affine Image Registration
            if app.RADIO_IR_LKFlowAffine.Value
                [LR_reg, Tvec, iter, err] = RegisterImageSeqAffine(app, stack);
            end
            
            %% LK Optical Flow Image Registration
            if app.RADIO_IR_LKFlow.Value
                [LR_reg, Tvec, iter, err] = RegisterImageSeq(app, stack);
            end
            
            % compute the mean number of iterations (rounded) and the
            % mean relative error
            iter = ceil(iter/(size(stack,3)-1));
            err = err/(size(stack,3)-1);
        end
        
        % Computes the HR frame from the registered LR stack
        function [HR, iter] = superResolution(app, LR, Tvec, resFactor, Hpsf, params)
            ShowProgress(app, ' Estimating High Resolution image', 0);
            
            %% Adaptive Kernel Regression
            if app.RADIO_SR_KernelReg.Value
                [HR, iter] = AdaptiveKernelRegression(app, LR, resFactor, params);
            end
            
            %% Cubic Spline Interpolation
            if app.RADIO_SR_Cubic.Value
                [HR, iter] = CubicSplineInterp(app, LR, resFactor, Hpsf, params);
            end
            
            %% Robust SR
            if app.RADIO_SR_Robust.Value
                [HR, iter] = RobustSR(app, LR, Tvec, resFactor, Hpsf, params);
            end
            
            %% Fast Robust SR
            if app.RADIO_SR_FastRobust.Value
                [HR, iter] = FastRobustSR(app, LR, Tvec, resFactor, Hpsf, params);
            end
            ShowProgress(app, '', 100);
        end

    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            
            % Load all components
            addpath([pwd '/Helper']);
            % Image Registration
            addpath([pwd '/ImageRegistration/Matlab']);
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
            
            % Set file formats allowed to open
            filter = '*.avi;*.mov;*.mp4;.m4v;';
            [FileName,PathName] = uigetfile(filter,'Select the movie file (avi, mov, mp4, m4v)');
            
            % Fix for GUI to get focus after loading file
            drawnow;
            figure(app.MFSRToolUIFigure);
            
            % Check if file exists
            if FileName ~= 0
                % Load video file
                inVideo = LoadVideo([PathName FileName]);
                
                % Set the LR stack
                app.LRstack = inVideo;
                app.frameCnt = 1;
                
                % Set preview image for LR stack
                app.IMG_LR.ImageSource = BuildImage(inVideo(:,:,1));
                
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
            app.IMG_LR.ImageSource = BuildImage(app.LRstack(:,:,app.frameCnt));
            
            % When passing the last frame, disable 'next' button
            if app.frameCnt == size(app.LRstack,3)
                app.BTN_next.Enable = 'off';
            end
            
            % When passing the first frame, enable 'previus' button
            if app.frameCnt > 1
                app.BTN_prev.Enable = 'on';
            end
        end

        % Button pushed function: BTN_prev
        function BTN_prevButtonPushed(app, event)
            
            % Display the previous frame in the stack
            app.frameCnt = app.frameCnt - 1;
            app.IMG_LR.ImageSource = BuildImage(app.LRstack(:,:,app.frameCnt));
            
            % At the first frame, disable 'previous' button
            if app.frameCnt == 1
                app.BTN_prev.Enable = 'off';
            end
            
            % If multiple frames, enable 'next' button
            if app.frameCnt < size(app.LRstack,3)
                app.BTN_next.Enable = 'on';
            end
        end

        % Button pushed function: BTN_saveLRFrame
        function BTN_saveLRFrameButtonPushed(app, event)

            % Save LR image to file
            image = app.LRstack(:,:,app.frameCnt);
            SaveFile(app, image);

        end

        % Button pushed function: BTN_enhance
        function BTN_enhanceButtonPushed(app, event)
            
            % Collect the parameters from the GUI
            [params, resFactor, Hpsf] = collectParams(app);
            
            % Compute the affine transformation matrices for each frame of
            % the sequence using the method chosen in the GUI
            if(~app.imReg_flag)
                tic;
                % Perform Image Registration
                [app.LR_reg, app.TM, iter, err] = getTransform(app);
                
                % Fill the benchmark measurement table
                app.VAL_IR_t.Text = num2str(toc);
                app.VAL_IR_err.Text = num2str(err);
                app.VAL_IR_n.Text = num2str(iter);
                
                % Set the image registration flag
                app.imReg_flag = true;
            end
            
            tic;
            % Perform Super Resolution
            [app.HRimage, iter] = superResolution(app, app.LRstack, app.TM, resFactor, Hpsf, params);
            
            % Fill the benchmark measurement table
            app.VAL_SR_t.Text = num2str(toc);
            app.VAL_SR_err.Text = 'n.A.'; % num2str(err);
            app.VAL_SR_n.Text = num2str(iter);
            
            % Set preview image of HR frame
            app.IMG_HR.ImageSource = BuildImage(app.HRimage);
            
            % Enable save button
            app.BTN_saveHRFrame.Enable = true;
            
        end

        % Button pushed function: BTN_saveHRFrame
        function BTN_saveHRFrameButtonPushed(app, event)
            
            % Save HR image to file
            image = app.HRimage;
            SaveFile(app, image);
            
        end

        % Selection changed function: RADIOGROUP_IR_Method
        function RADIOGROUP_IR_MethodSelectionChanged(app, event)
            
            % Reset image registration flag
            app.imReg_flag = false;
            
        end

        % Button pushed function: BTN_reset
        function BTN_resetButtonPushed(app, event)
            
            % Reset image registration flag
            app.imReg_flag = false;
            
            % Reset preview image of HR frame
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
            app.LABEL_imInfo_imHeight.FontSize = 13;
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
            app.LABEL_imInfo_imWidth.FontSize = 13;
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
            app.LABEL_SR_time.Text = 'time (sec):';

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
            app.LABEL_SR_err.Text = 'error:';

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

            % Create ResFactorLabel
            app.ResFactorLabel = uilabel(app.PANEL_SR_params);
            app.ResFactorLabel.FontSize = 14;
            app.ResFactorLabel.FontColor = [0.9412 0.9412 0.9412];
            app.ResFactorLabel.Position = [10 57 90 22];
            app.ResFactorLabel.Text = 'Res-Factor';

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