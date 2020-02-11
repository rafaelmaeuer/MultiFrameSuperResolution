function SaveFile(app, file)
    
    % Set file formats allowed to save
    [FileName,PathName] = uiputfile('*.jpg','Save image file');

    % Fix for GUI to get focus after loading file
    drawnow;
    figure(app.MFSRToolUIFigure);

    % If filename exists save file
    if FileName ~= 0
      imwrite(uint8(file), [PathName FileName]);
    end
    
end

