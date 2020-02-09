function [croppedIm] = APP_cropImage(theImage, refWidth, refHeight)
    
    % get dimensions of the image to be cropped
    imHeight = size(theImage,1);
    imWidth = size(theImage,2);

    croppedIm (1:refHeight, 1:refWidth) = 128;
    
    if imHeight < refHeight
        y1 = floor((refHeight-imHeight)/2)+1;
        y2 = 1;
        h = imHeight -1;
        
    else
        y1 = 1;
        y2 = floor((imHeight-refHeight)/2)+1;
        h = refHeight -1;
    end
    
    if imWidth < refWidth
        x1 = floor((refWidth-imWidth)/2)+1;
        x2 = 1;
        w = imWidth - 1;
        
    else
        x1 = 1;
        x2 = floor((imWidth-refWidth)/2)+1;
        w = refWidth - 1;
    end
    
    croppedIm (y1:y1+h, x1:x1+w) = theImage(y2:y2+h,x2:x2+w);
end

