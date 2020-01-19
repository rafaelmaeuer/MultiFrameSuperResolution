% Loads the given AVI file and convert it to an image struct
function M=LoadVideo(filename)

videoIn= VideoReader(filename);

s = 1;

while hasFrame(videoIn)
    
    tempBuffer = readFrame(videoIn);
    buffer(:,:,s)= double(tempBuffer(:,:,1));
    s = s+1;
end

M=buffer;

if (length(size(M))==4)
  M=squeeze(M(:,:,1,:));
end