% Loads the given MAT file which contains a video file
% If the video is a color image, then it is converted to gray
function M=LoadVideoMat(filename)

s=load(filename);
s=struct2cell(s);
M=double(s{1});

if (length(size(M))==4)
  M=squeeze(M(:,:,1,:));
end