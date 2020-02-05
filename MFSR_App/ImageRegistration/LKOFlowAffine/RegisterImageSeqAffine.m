% D=RegisterImageSeqAffine(M)
%
% Affine Image-Registration of a given sequence of images
function D=RegisterImageSeqAffine(M)

% Initialize d to an empty affine displacement transformation
D=zeros(2,3,size(M,3));

roi=[2 2 size(M,1)-1 size(M,2)-1];

Mprev = squeeze(M(:,:,1));

h=waitbar(0, 'Calculating Affine Registration');

d=[1 0 0; 0 1 0];
D(:,:,1)=d;

for i=2:size(M,3)
  
  waitbar(i/size(M,3));
  
  % Register current image to previous frame
  dc=PyramidalLKOpticalFlowAffine(Mprev, squeeze(M(:,:,i)), roi);
  
  % Save current image
  Mprev = squeeze(M(:,:,i));

  % Add current displacement to d (This is actually concatinating the two
  % affine matrixes)
  d=d+reshape(dc, 2, 3)*(eye(3)+[d;0 0 0]);

  % Compute displacement at current level
  d=IterativeLKOpticalFlowAffine(squeeze(M(:,:,1)), squeeze(M(:,:,i)), roi, d);
  
  % Transform current image
  %M(:,:,i)=ResampleImgAffine(M(:,:,i), [1 1 size(M,1) size(M,2)], d);
  D(:,:,i)=d;

  %figure;imagesc(M(:,:,i));title('reg2base refined');set(gcf,'name', 'reg2base refined');

end

close(h);