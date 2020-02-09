classdef affine_flow
    % AFFINE_FLOW estimates affine (first-order) optic flow.
    %
    % Given two images, this class computes estimates of the parameters of
    % an affine (first-order) flow model which describes the overall
    % deformation between the two images.
    %
    % Affine flow
    % -----------
    %
    % An affine flow field has six parameters. One way to express these is
    % as follows. If the optic flow vector at an image location (x,y) is
    % (vx,vy), the first-order model is:
    %
    %           [vx vy] = [x y 1] [  d+s1  s2+r  ]
    %                             [  s2-r  d-s1  ]
    %                             [   vx0   vy0  ]
    %
    % Note that when x and y refer to pixels in the image array, the
    % conventional order is reversed, and the image intensity at (x,y) is
    % given by IMAGE(y,x) - that is, y is the row and x is the column.
    %
    % The parameters may be thought of thus: (vx0, vy0) is the optic flow
    % at the origin (one pixel outside the top left corner of the image);
    % d is the rate of dilation; r is the rate of rotation; s1 is the
    % shear along the main image axes; s2 is the shear along the diagonal
    % axes.
    %
    % Algorithm
    % ---------
    %
    % The flow is estimated by carrying out a least-squares fit of the
    % model to the spatial and temporal gradients of the image pair; this
    % is an extension of the Lucas-Kanade method. Gradients are estimated
    % with simple symmetrical differences and optional smoothing.
    %
    % Optionally, the images may be resampled on a log-polar grid, as
    % described here: http://www.bmva.org/bmvc/1994/bmvc-94-057.pdf. This
    % may increase the estimation accuracy because the grid spacing can be
    % matched to the likely magnitude of the optic flow. When log-polar
    % sampling is used, an iterative scheme to refine the estimates by
    % moving the sampling centre with the flow may be activated.
    %
    % The simplest use of the class, with no properties set explicitly, is
    % unlikely to yield good results on real images. At a minimum,
    % smoothing of the images is needed during computation of the
    % gradients. For more information, see the help for the individual
    % properties.
    %
    % Usage
    % -----
    %
    % This is a value class, so any update to the state of an affine_flow
    % object requires an assignment.
    %
    % A new affine_flow object is set up using parameter/value arguments,
    % for example:
    %
    %   af = affine_flow('sigmaXY', 10, 'sampleStep', 2);
    %
    % Parameters may also be updated subsequently. Parameters that are not
    % set are given defaults. To see the property names and default values
    % for rectilinear sampling, do
    %
    %   disp(affine_flow);
    %
    % To see the property names and default values for log-polar sampling,
    % do
    %
    %   disp(affine_flow('sampleMethod', 'logpolar'));
    %
    % The images to be processed may be supplied to the constructor, or
    % subsequently:
    %
    %   af.image1 = im1;
    %   af.image2 = im2;
    %
    % The flow is computed by a call to the findFlow method:
    %
    %   af = af.findFlow;
    %
    % and the results obtained from the flowStruct property:
    %
    %   flow = af.flowStruct;
    %
    % See the individual properties and methods for more details, and
    % information about the results structure.
    %
    % Efficiency notes
    % ----------------
    %
    % If many pairs of images are to be processed, it is most efficient to
    % create a single affine_flow object outside the loop and use it to
    % process each pair of images.
    %
    % If multiple images are to be matched to a given initial image, only
    % the image2 property should be updated within the loop, to avoid
    % refiltering the first image.
    %
    % If the flow between successive pairs of a sequence of images is
    % required, it is most efficient to call the advance method within the
    % loop to move image2 to image1, and only to update the image2 property
    % directly.
    %
    % See also: AFFINE_FLOWDEMO, LOGSAMPLE, GRADIENTS_XYT, 
    % AFFINE_FLOWDISPLAY, AFFINE_FLOWEDGEDISPLAY

    % Copyright David Young 2010

    
    properties (Dependent)
        
        % Sampling method
        %
        % Values are 'rectilinear' or 'logpolar'; the default is
        % 'rectilinear'. The image gradients are sampled either on a
        % regular rectilinear grid aligned with the image axes, or on a
        % log-polar grid. The possible advantage of a log-polar grid is
        % that the sample spacing increases with distance from the centre,
        % so provided the optic flow is small at the grid centre, or the
        % grid is moved to track the motion, the sample spacing will be
        % better matched to the flow velocity variation across the image.
        %
        % See also logsample
        sampleMethod
        
        % General region of interest
        %
        % A binary array (double or logical) the same size as the images,
        % used as a mask, or the empty matrix. Gradients are only computed
        % for pixels where there is a non-zero value in the roi array. The
        % smoothing process may access pixels outside the ROI mask. Ignored
        % if set to the empty matrix (the default).
        %
        % There are simpler ways to restrict processing to part of the
        % image: for rectilinear sampling, a rectangular region may be
        % specified with the rectRegion property; for log-polar sampling, a
        % circular region may be specified with logCentre and logRmax.
        %
        % Setting regionOfInterest will discard any previous value for the
        % rectRegion property. These two properties should not both be set
        % in a call to the class constructor.
        regionOfInterest
        
        % Rectangular region of interest (rectilinear sampling only)
        %
        % A 4-element row vector of integers, or the empty matrix. This
        % specifies a rectangular region of the image in which to compute
        % the flow, in the form [ROWMIN ROWMAX COLMIN COLMAX]. Because this
        % resembles array indexing, (row,column) rather than (x,y) order is
        % used. Gradients are computed only for pixels in the region, but
        % the smoothing process may access pixels just outside the region. 
        % 
        % If set to the empty matrix (the default), the region used is that
        % for which the smoothing and differencing convolutions, applied to
        % the whole image, produce 'valid' results (i.e. results which do
        % not require a guess at what lies outside the image boundary).
        %
        % Setting rectRegion will discard any previous value for the
        % regionOfInterest property. These two properties should not both
        % be set in a call to the class constructor.
        rectRegion
        
        % Smoothing constant for rectilinear sampling
        %
        % A non-negative real scalar; default 0. The parameter for the
        % Gaussian mask used to smooth the images for gradient estimation.
        %
        % See also GSMOOTH2
        sigmaXY
        
        % Sampling interval (rectilinear sampling only)
        %
        % A positive integer; default 1. The system of equations that is
        % solved is heavily overdetermined and it is not generally
        % necessary to use gradients for every pixel for the least-squares
        % solution. Gradients can therefore be sampled on a grid with a
        % spacing of sampleStep pixels in x and y. The sampling is done
        % after smoothing. It is usually reasonable to set sampleStep
        % roughly equal to sigmaXY.
        sampleStep

        % Smoothing constants for log-polar sampling
        %
        % Either a non-negative real scalar or a 2-element row vector of
        % non-negative real numbers; default 0. If a vector, the first
        % element is the parameter for the Gaussian mask used to smooth the
        % log-polar sampled image across rings (i.e. for a mask oriented
        % along the ring axis); the second element is the parameter for
        % smoothing across wedges. A scalar sets both parameters.
        sigmaRW
     
        % Centre of log-polar sampling pattern
        %
        % A 2-element positive real row vector of the form [X0 Y0], or the
        % empty matrix. X0 specifies the x (column) coordinate of the
        % centre of the log-polar sampling grid. Y0 specifies the y (row)
        % coordinate of the centre. If empty (the default) the centre of
        % the image is used.
        %
        % See also logsample
        logCentre
        
        % Inner radius of log-polar sampling pattern
        %
        % A positive integer; default 1. The radius of the smallest ring of
        % the log-polar sampling grid.
        %
        % See also logsample
        logRmin

        % Outer radius of log-polar sampling pattern
        %
        % A positive integer or the empty matrix, specifying the radius of
        % the largest ring of the log-polar sampling grid. If empty (the
        % default) and maxIter is 1, the largest value that keeps the outer
        % ring entirely within the image is used; if maxIter is greater
        % than 1, defaults to 75% of this.
        %
        % See also logsample
        logRmax

        % Number of wedges in the log-polar sampling pattern
        %
        % A positive integer or the empty matrix, specifying the number of
        % wedges in the log-polar sampling grid. If empty (the default)
        % round(logRmax) is used; this gives an arbitrary sample spacing of
        % 2*pi pixels in the outer ring.
        %
        % The number of rings is chosen automatically to satisfy the
        % "circular samples" constraint.
        %
        % See also logsample
        logWedges

        % Maximum number of iterations (log-polar sampling only)
        %
        % A positive integer; default 1. If greater than 1 the log-polar
        % grid applied to image2 is moved to track the current estimated
        % optic flow for the grid centre, and the flow is then reestimated.
        % This can effectively null out the flow at the centre of the grid,
        % allowing more accurate estimation. The flow parameters returned
        % are corrected for the grid motion.
        %
        % Iteration stops when the number of iterations reaches maxIter, or
        % when the residual velocity drops below velTol, or when the outer
        % ring of the log-polar grid would go outside the image.
        %
        % If regionOfInterest has been specified, the ROI for the second
        % image moves with the log-polar grid during tracking.
        maxIter

        % Speed estimation tolerance for tracking (log-polar sampling only)
        %
        % A positive real scalar; default 1. If tracking is in effect,
        % iteration will halt when the magnitude of the residual flow at
        % the centre of the grid has dropped below velTol. The residual
        % flow is the estimated true flow minus the offset of the log-polar
        % grid between IM1 and IM2.
          velTol
        
        % The initial image
        %
        % Must be set by assignment or as an argument to the constructor
        % before findFlow can be called. Must be set to a 2D array.
        image1
        
        % The image after motion has taken place
        %
        % Must be set by assignment or as an argument to the constructor
        % before findFlow can be called. Must be set to a 2D array. Must
        % have the same size as image1.
        image2
        
    end
    
    properties (Dependent, SetAccess = private)
        
        % The affine flow estimate as a structure
        %
        % The estimated parameters. The structure has fields vx0, vy0, d,
        % r, s1 and s2. For the meanings of the parameters, see the
        % introductory section of the help for the class.
        flowStruct
        
        % The affine flow estimate as a matrix
        %
        % The estimated parameters, stored as a matrix M such that the
        % estimated flow at a position (x, y) is given by [x y 1]*M.
        flowMatrix
        
        % The estimate of the image deformation matrix
        %
        % The estimated parameters, stored as a matrix W such that if (x,
        % y) is a position in image1, the estimate of the corresponding
        % position in image2 is given by [x y 1] * W.
        flowWarp
        
        % Tracking information (log-polar sampling only)
        %
        % The effect of tracking during the last call to findFlow is
        % reported in the trackInfo structure, whose fields are:
        %
        %   numIter - the number of iterations that took place;
        %   lastCentre - a 2-element row vector containing the x and y
        %       coordinates of the centre of the final sampling pattern
        %       applied to image2;
        %   lastspeed, the final residual flow speed at the log-polar
        %       centre after tracking has been applied.
        trackInfo
        
    end
    
    properties (Access = private)
        
        % Shadows of dependent variables
        
        im1 = []
        im2 = []
        rectmeth = true
        roi = []
        region = []
        sigmaxy = 0
        sampstep = 1
        sigmarw = [0 0]
        x0 = []
        y0 = []
        rmin = 1
        rmax = []
        nw = []
        maxiter = 1
        veltol = 1
        
        % Intermediate results, held to possibly reduce later computation
        
        smth1 = []
        smth2 = []
        roireg = []
        roitrimmed = []
        
        logsmth1 = []
        logsmth2 = []
        logreg = []
        logroi = []
        
        % Results of main computation
        
        flow = []
        resOK = false
        iterations = []
        lastx0 = []
        lasty0 = []
        lastspeed = []
    end
    
    
    methods         % constructor and get/set methods
        
        function a = affine_flow(varargin)
            % Constructor for affine flow objects
            %
            % af = affine_flow('prop1', val1, 'prop2', val2, ...)
            % constructs an affine flow object, setting property prop1 to
            % value val1 etc. See the class information for details.
            % Properties not set are given default values; to see these
            % execute disp(affine_flow) and
            % disp(affine_flow('sampleMethod', 'logpolar')).
            %
            % See also: affine_flowdemo
            a = setProps(a, varargin{:});
        end
        
        % Set/get methods with argument checking and state maintenance.
        % Care is needed to ensure that intermediate results are correctly
        % (but not unnecessarily) invalidated if properties are changed.
        
        function a = set.image1(a, im)
            % Could validate image1, but don't for speed
            a.im1 = im;
            a.smth1 = [];
            a.logsmth1 = [];
            a.resOK = false;
        end
        function im = get.image1(a)
            im = a.im1;
        end
        
        function a = set.image2(a, im)
            % Could validate image2, but don't for speed
            a.im2 = im;
            a.smth2 = [];
            a.logsmth2 = [];
            a.resOK = false;
        end
        function im = get.image2(a)
            im = a.im2;
        end
        
        function a = set.sampleMethod(a, s)
            if ~ismember(s, {'rectilinear' 'logpolar'})
                error('affine_flow:badmethod', ...
                    'Method must be rectilinear or logpolar');
            end
            a.rectmeth = isequal(s, 'rectilinear');
            a.resOK = false;
        end
        function s = get.sampleMethod(a)
            if a.rectmeth
                s = 'rectilinear';
            else
                s = 'logpolar';
            end
        end
        
        function a = set.regionOfInterest(a, roi)
            if ~isempty(roi)
                validateattributes(roi, {'double' 'logical'}, {'binary'});
            end
            a.roi = roi; a.roireg = []; a.logroi = [];
            % rect preprocessing needs to be redone because bounding box
            % may have changed, but lp preprocessing doesn't
            a.smth1 = []; a.smth2 = [];
            a.region = [];
            a.resOK = false;
        end
        function r = get.regionOfInterest(a)
            r = a.roi;
        end
        
        function a = set.rectRegion(a, reg)
            if ~isempty(reg)
                validateattributes(reg, {'double'}, {'integer' 'size' [1 4]});
            end
            a.roi = []; a.roireg = []; a.logroi = [];
            a.smth1 = []; a.smth2 = [];
            a.region = reg;
            a.resOK = false;
        end
        function r = get.rectRegion(a)
            a.checkrelevance('rectRegion', true);
            r = a.region;
        end
        
        function a = set.sigmaXY(a, sigma)
            validateattributes(sigma, {'double'}, ...
                {'nonnegative' 'real' 'finite' 'scalar'});
            a.sigmaxy = sigma;
            a.smth1 = []; a.smth2 = [];
            a.resOK = false;
        end
        function s = get.sigmaXY(a)
            a.checkrelevance('sigmaXY', true);
            s = a.sigmaxy;
        end
        
        function a = set.sampleStep(a, step)
            validateattributes(step, {'double'}, ...
                {'positive' 'integer' 'scalar'});
            a.sampstep = step;
            a.resOK = false;
        end
        function s = get.sampleStep(a)
            a.checkrelevance('sampleStep', true);
            s = a.sampstep;
        end
        
        function a = set.sigmaRW(a, sigma)
            validateattributes(sigma, {'double'}, ...
                {'nonnegative' 'real' 'finite'});
            if isscalar(sigma)
                a.sigmarw = [sigma sigma];
            else
                if ~isequal(size(sigma), [1 2])
                    error('affine_flow:badsigmaRW', ...
                        'sigmaRW must be scalar or 1x2 vector');
                end
                a.sigmarw = sigma;
            end
            a.logsmth1 = []; a.logsmth2 = []; a.logroi = [];
            a.resOK = false;
        end
        function s = get.sigmaRW(a)
            a.checkrelevance('sigmaRW', false);
            s = a.sigmarw;
        end
        
        function a = set.logCentre(a, c)
            if isempty(c)
                a.x0 = [];
                a.y0 = [];
            else
                validateattributes(c, {'double'}, ...
                    {'real' 'finite' 'scalar' 'size' [1 2]});
                a.x0 = c(1);
                a.y0 = c(2);
            end
            a.logsmth1 = []; a.logsmth2 = []; a.logroi = [];
            a.resOK = false;
        end
        function c = get.logCentre(a)
            a.checkrelevance('logCentre', false);
            c = [a.x0 a.y0];
        end
        
        function a = set.logRmin(a, r)
            validateattributes(r, {'double'}, ...
                {'positive' 'real' 'finite' 'scalar'});
            a.rmin = r;
            a.logsmth1 = []; a.logsmth2 = []; a.logroi = [];
            a.resOK = false;
        end
        function r = get.logRmin(a)
            a.checkrelevance('logRmin', false);
            r = a.rmin;
        end
        
        function a = set.logRmax(a, r)
            if ~isempty(r)
                validateattributes(r, {'double'}, ...
                    {'positive' 'real' 'finite' 'scalar'});
            end
            a.rmax = r;
            a.logsmth1 = []; a.logsmth2 = []; a.logroi = [];
            a.resOK = false;
        end
        function r = get.logRmax(a)
            a.checkrelevance('logRmax', false);
            r = a.rmax;
        end
        
        function a = set.logWedges(a, n)
            validateattributes(n, {'double'}, {'positive' 'integer' 'scalar'});
            a.nw = n;
            a.logsmth1 = []; a.logsmth2 = []; a.logroi = [];
            a.resOK = false;
        end
        function n = get.logWedges(a)
            a.checkrelevance('logWedges', false);
            n = a.nw;
        end
        
        function a = set.maxIter(a, n)
            validateattributes(n, {'double'}, {'positive' 'integer' 'scalar'});
            a.maxiter = n;
            a.resOK = false;
        end
        function n = get.maxIter(a)
            a.checkrelevance('maxIter', false);
            n = a.maxiter;
        end
        
        function a = set.velTol(a, v)
            validateattributes(v, {'double'}, ...
                {'positive' 'real' 'finite' 'scalar'});
            a.veltol = v;
            a.resOK = false;
        end
        function v = get.velTol(a)
            a.checkrelevance('velTol', false);
            v = a.veltol;
        end
        
        function f = get.flowStruct(a)
            if a.resOK
                f = a.flow;
            else
                error('affine_flow:noresult', 'No flow result available');
            end
        end
        
        function m = get.flowMatrix(a)
            m = affine_flow.matrix(a.flowStruct);
        end
        
        function w = get.flowWarp(a)
            w = affine_flow.warp(a.flowStruct);
        end
        
        function t = get.trackInfo(a)
            if a.resOK && ~a.rectmeth
                t.lastCentre = [a.lastx0 a.lasty0];
                t.lastSpeed = a.lastspeed;
                t.numIter = a.iterations;
            else
                error('affine_flow:notrack', 'No tracking result available');
            end
        end
            
    end
    
    
    methods         % public computational methods
        
        function a = findFlow(a)
            % Computes affine flow between image1 and image2
            %
            % Estimates flow using extended Lucas-Kanade method (see
            % algorithm description). Results are accessed via the
            % flowStruct/flowMatrix/flowWarp and trackInfo properties.
            
            if a.rectmeth
                a = a.affine_flowxy;
            else
                a = a.affine_flowrw;
            end
        end
        
        function a = advance(a)
            % advance advances image2 to image1.
            %
            % This method is for computation of successive sets of
            % flow fields in image sequences, by moving image2 to
            % image1 and updating the feature positions to the new
            % positions in image2. Its main effect is thus
            %
            %   af.image1 = af.image2;
            %
            % However, the method is more efficient than the assignment
            % would be. It could be used like this, assuming we have a cell
            % array of images and we want a cell array of flow field
            % results:
            %
            %   af = affine_flow('image1', images{1});
            %   for k = 2:length(images)
            %       af.image2 = images{k};
            %       af = af.findFlow;
            %       flowfields{i} = af.flowStruct;
            %       af = af.advance;
            %   end
            %
            % If log-polar sampling is used and tracking is in effect, the
            % initial sampling centre also moves to the final centre after
            % tracking, in effect:
            %
            %   af.logCentre = af.trackInfo.lastCentre;
            
            a.im1 = a.im2;
            a.smth1 = a.smth2;
            a.logsmth1 = a.logsmth2;
            a.x0 = a.lastx0;
            a.y0 = a.lasty0;
            a.resOK = false;
        end
        
    end
    
    
    methods (Access = private)
        
        function a = affine_flowxy(a)
            % Main computational function for rectilinear sampling.
            
            roil = a.roi;
            regionl = a.region;
            
            if ~isempty(roil)
                if isempty(a.roireg)
                    % Trim all arrays to the bounding box to avoid
                    % unnecessary work. Convert to double for regionprops
                    % as region may be unconnected
                    bbox = regionprops(double(a.roi), 'BoundingBox');
                    bbox = bbox.BoundingBox;
                    regionl = [ceil(bbox(2)) floor(bbox(2)+bbox(4)) ...
                        ceil(bbox(1)) floor(bbox(1)+bbox(3))];
                    roil = logical(a.getreg(a.roi, regionl));
                    a.roireg = regionl;
                    a.roitrimmed = roil;
                else
                    regionl = a.roireg;
                    roil = a.roitrimmed;
                end
            elseif isempty(regionl)
                regionl = gradients_xyt(a.im1, a.im2, a.sigmaxy, ...
                    'region', false);
            end
            
            % Image gradients and x,y coordinate arrays
            if isempty(a.smth1)
                a.smth1 = gsmooth2(a.im1, a.sigmaxy, ...
                    regionl + [-1 1 -1 1], false);
            end
            if isempty(a.smth2)
                a.smth2 = gsmooth2(a.im2, a.sigmaxy, ...
                    regionl + [-1 1 -1 1], false);
            end
            [xg, yg, tg] = gradients_xyt(a.smth1, a.smth2, 0, [], false);
            
            % Use origin at centre of image for now - more stable
            [x, y, cx, cy] = a.centredxy(xg);
            
            % Cut down typing by assembling the arguments
            flowargs = {x, y, xg, yg, tg};
            
            % Resample on bigger grid if required to cut down least-squares
            % solver effort
            if a.sampstep > 1
                steps = [a.sampstep a.sampstep];
                flowargs = cellfun(@(z) {affine_flow.sample(z, steps)}, ...
                    flowargs);
            end
            
            % Get values in the ROI if there is one, otherwise just convert
            % to column vectors
            if ~isempty(roil)
                if a.sampstep > 1
                    roil = affine_flow.sample(roil, steps);
                end
                flowargs = cellfun(@(z) {z(roil)}, flowargs);
            else
                flowargs = cellfun(@(z) {z(:)}, flowargs);
            end
            
            % Assemble and solve the least squares system
            flowl = affine_flow.solvexy(flowargs{:});
            
            % Move origin back to origin of image coordinates
            a.flow = affine_flow.shift(flowl, ...
                1-regionl(3)-cx, 1-regionl(1)-cy);
            a.resOK = true;
            
        end
        
        
        function a = affine_flowrw(a)
            % Main computational function for log-polar sampling.
            
            centre = (size(a.im1)+1)/2;
            xc = a.x0;
            if isempty(xc)
                xc = centre(2);
            end
            yc = a.y0;
            if isempty(yc)
                yc = centre(1);
            end
            rmx = a.rmax;
            if isempty(rmx)
                % keep samples just within image
                rmx = min([centre-1 size(a.im1)-centre]);
                % but if we are iterating, allow room for movement
                if a.maxIter > 1
                    rmx = 0.75 * rmx;   % arbitrary
                end
            end
            nwed = a.nw;
            if isempty(nwed)
                nwed = round(rmx);     % quite arbitrary
            end
            
            [maxy, maxx] = size(a.im1);
            
            % Sample first image onto log-polar grid
            if isempty(a.logsmth1)
                logim1 = logsample(a.im1, a.rmin, rmx, xc, yc, [], nwed);
                [a.logsmth1, regions] = gsmooth2(logim1, a.sigmarw, [], ...
                    [false true]);
                a.logreg = regions + [0 0 1 -1];
            end
            
            [r, w] = affine_flow.logcoords(a.logreg);
            
            % and roi
            roil = a.roi;
            if isempty(a.logroi) && ~isempty(roil)
                % imtransform (and so logsample) is OK with logical images
                logroil = logsample(logical(roil), a.rmin, rmx, ...
                    xc, yc, [], nwed);
                a.logroi = affine_flow.getreg(logroil, a.logreg);
            end
            
            xc2 = xc;
            yc2 = yc;
            iter = 0;
            
            while true      % iterate to track motion
                iter = iter + 1;
                
                % Sample second image onto current centre for second grid
                if isempty(a.logsmth2)
                    logim2 = logsample(a.im2, a.rmin, rmx, xc2, yc2, [], nwed);
                    a.logsmth2 = gsmooth2(logim2, a.sigmarw, [], ...
                        [false true]);
                end
                
                % Get gradients
                [gr, gw, gt] = gradients_xyt(a.logsmth1, a.logsmth2, ...
                    0, [], [false true]);
                
                % Cut down typing by assembling the arguments
                flowargs = {r, w, gr, gw, gt};
                
                % Get values in the ROI if there is one, otherwise just
                % convert to column vectors
                if ~isempty(roil)
                    flowargs = cellfun(@(z) {z(a.logroi)}, flowargs);
                else
                    flowargs = cellfun(@(z) {z(:)}, flowargs);
                end
                
                % Build and solve the least-squares system
                flowl = affine_flow.solverw(a.rmin, nwed, flowargs{:});
                
                % Prepare for tracking if required. save most recent sample
                % centre
                a.lastx0 = xc2;
                a.lasty0 = yc2;
                % update centre position to follow flow
                xc2 = xc2 + flowl.vx0;
                yc2 = yc2 + flowl.vy0;
                % test whether finished
                speed = norm([flowl.vx0 flowl.vy0]);
                if iter >= a.maxiter || speed <= a.veltol || ...
                        min([xc2-1 yc2-1 maxx-xc2 maxy-yc2]) < rmx
                    break;
                end
                
                a.logsmth2 = [];
            end
            
            % Correct centre velocity for sample centre shift
            flowl.vx0 = xc2 - xc;
            flowl.vy0 = yc2 - yc;
            % Move origin back to origin of image coordinates
            a.flow = affine_flow.shift(flowl, -xc, -yc);
            a.resOK = true;
            
            % Record extra information
            a.iterations = iter;
            a.lastspeed = speed;
            
        end
        
        function checkrelevance(a, propname, rectprop)
            if rectprop ~= a.rectmeth
                error('affine_flow:badfieldrequest', ...
                    ['Property ' propname ...
                    ' not used for current sampling method']);
            end
        end
        
    end
    
    
    methods (Static, Access = private)
        
        function [x, y, cx, cy] = centredxy(im)
            % Makes x and y arrays centred on the centre of the image
            [s1, s2] = size(im);
            cy = (s1-1)/2;
            cx = (s2-1)/2;
            [x, y] = meshgrid(linspace(-cx, cx, s2), linspace(-cy, cy, s1));
        end
        
        
        function [r, w] = logcoords(reg)
            % makes r and w index arrays for the specified region of the
            % log-sampled image
            [r, w] = meshgrid(reg(3)-1:reg(4)-1, reg(1)-1:reg(2)-1);
        end
        
        
        function arr = getreg(arr, reg)
            % trim array to specified region
            arr = arr(reg(1):reg(2), reg(3):reg(4));
        end
        
        
        function arr = sample(arr, steps)
            % sample every sx'th column and sy'th row, centering the grid
            sz = size(arr);
            nm1 = floor((sz-1)./steps);
            starts = floor((sz - (steps.*nm1 + 1))/2) + 1;
            arr = arr(starts(2):steps(2):end, starts(1):steps(1):end);
        end
        
        
        function affineparams = solvexy(x, y, gx, gy, gt)
            % AFFINE.SOLVEXY first order optic flow parameters -rectilinear
            %   AFFINEPARAMS = AFFINE.SOLVEXY(X, Y, GXS, GYS, GTS) returns
            %   the affine flow parameters vx0, vy0, d, r, s1, s2, as a
            %   struct with those f ields. The inputs are the x and y
            %   sample positions and the x, y and t gradients, all as
            %   column vectors.
            
            xgx = x .* gx;
            ygx = y .* gx;
            xgy = x .* gy;
            ygy = y .* gy;
            A = [gx, gy, xgx+ygy, xgy-ygx, xgx-ygy, ygx+xgy];
            
            a = -(A \ gt);
            
            affineparams = struct('vx0', a(1), 'vy0', a(2), ...
                'd', a(3), 'r', a(4), 's1', a(5), 's2', a(6));
        end
        
        
        function affineparams = solverw(rmin, nw, r, w, gr, gw, gt)
            % AFFINE.SOLVERW first order optic flow parameters - log-polar
            %   AFFINEPARAMS = AFFINE.SOLVERW(RMIN, NW, R, W, GR, GW, GT)
            %   returns the affine flow parameters vx0, vy0, d, r, s1, s2,
            %   as a struct with those fields. The inputs are the
            %   parameters of the log-polar grid and the gradients measured
            %   on it.
            
            kinv = 2*pi/nw;
            rhoinv = exp(-kinv*r)/rmin;
            theta = w * kinv;
            c = cos(theta);
            s = sin(theta);
            grc = gr .* c;
            grs = gr .* s;
            gwc = gw .* c;
            gws = gw .* s;
            
            gx = (grc - gws) .* rhoinv;
            gy = (grs + gwc) .* rhoinv;
            gs1 = grc.*c - grs.*s - 2*gwc.*s;
            gs2 = 2*grc.*s + gwc.*c - gws.*s;
            
            A = [gx, gy, gr, gw, gs1, gs2];
            
            a = -kinv * (A \ gt);
            
            affineparams = struct('vx0', a(1), 'vy0', a(2), ...
                'd', a(3), 'r', a(4), 's1', a(5), 's2', a(6));
        end
        
    end
    
    
    methods (Static)
        
        function m = matrix(f)
            %AFFINE_FLOW.MATRIX converts affine flow structure to matrix
            %   M = AFFINE_FLOW.MATRIX(F) takes an affine flow structure as
            %   returned by AFFINE_FLOW and returns a matrix M such that if
            %   P is a row vector with components [X, Y, 1] representing a
            %   position, P*M is a vector representing the optic flow at
            %   that position.
            m = [f.d+f.s1,  f.s2+f.r;
                f.s2-f.r,  f.d-f.s1;
                f.vx0,     f.vy0];
        end
        
        function w = warp(f)
            %AFFINE_FLOW.WARP converts a flow to a warp
            %   W = AFFINE_FLOW.WARP(F) takes either a flow structure as
            %   returned by AFFINE_FLOW or a flow matrix as returned by
            %   AFFINE_FLOW2MATRIX and returns a matrix W such that if P is
            %   a row vector with components [X, Y, 1] representing a
            %   position, then P*M is a vector representing the new
            %   position of the vector after one frame of the flow has
            %   occured.
            %
            %   Uses a very simple approximation!
            if isstruct(f)
                f = affine_flow.matrix(f);
            end
            w = f + [eye(2); 0 0];
        end
        
        function f = shift(f, x0, y0)
            %AFFINE_FLOW.SHIFT shifts the origin of affine flow
            %   F = AFFINE_FLOW.SHIFT(F, X0, Y0) returns affine flow
            %   parameters as returned by AFFINE_FLOW, with the origin
            %   shifted to X0, Y0 relative to the current origin.
            f.vx0 = f.vx0 + x0*(f.d+f.s1) + y0*(f.s2-f.r);
            f.vy0 = f.vy0 + x0*(f.s2+f.r) + y0*(f.d-f.s1);
        end
        
    end
    
end
