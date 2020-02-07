function [croppedIm] = APP_cropImage(theImage, rect)
    croppedIm = theImage(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
end

