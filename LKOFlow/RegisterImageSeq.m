function D=RegisterImageSeq(M)

% Initialize d to an empty affine displacement transformation
D=zeros(size(M,3),2);

roi=[2 2 size(M,1)-1 size(M,2)-1];

Mbase = squeeze(M(:,:,1));

h=waitbar(0, 'Calculating Registration');

for i=2:size(M,3)
  
  waitbar(i/size(M,3));
  
  % Register current image to previous frame
  D(i,:)=PyramidalLKOpticalFlow(Mbase, squeeze(M(:,:,i)), roi);
  
end

close(h);