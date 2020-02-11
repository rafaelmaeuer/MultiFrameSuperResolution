function img = BuildImage(theFrame)
    % display first frame in LR Image container.
    % IMG object expects a MxNx3 (rgb) matrix, so we just copy the
    % pixel data to all color channels:

    height = size(theFrame,1);
    width = size(theFrame,2);

    img = zeros(height,width,3);

    for i=1:3
      img (:,:,i) = double(theFrame(:,:))./255;
    end
end

