function segment_virions(fname);

load([fname,'.mat']);
K = 5;
w = 1;
for k = 1:K
  % Convert image matrices into doubles, so that we can perform mathmatical
  % operations on them:
  chl555 = double(ABS(:,:,k)); chl488 = double(LABEL(:,:,k));
  %disp(chl488)
  % Create binary images from the HA and NA channels, and combine:
  chl555_0 = chl555 > (median(chl555(:))+0.3*std(chl555(:)));
  chl488_0 = chl488 > (median(chl488(:))+std(chl488(:))); 
  %T = adaptthresh(f,0.0);
  %f0 = imbinarize(f,T);
  %f0 = a0; f0 = f0>0;
  sum_mask = chl555_0+0*chl488_0; sum_mask = sum_mask>0;
  %imshow(uint16([100*chl488 2^16*sum_mask])); drawnow; %pause;
  %disp(sum(sum(sum_mask)))
  % Erode and dilate the binary image to remove stray / spurious pixels:
  se = strel('disk',1);
  f1 = imerode(sum_mask,se);
  m0 = imdilate(f1,se);
  m1 = bwlabel(m0);
  %imshow(uint16([100*chl555 2^16*m0])); drawnow; %pause;

  % Using REGIONPROPS, extract information for each object (virus particle)
  % in the binary image. Crop the image around each separate object using
  % the bounding box.
  stats = regionprops(m1,"BoundingBox","MajorAxisLength","MinorAxisLength");
  size(stats);
  %MajorAxis=zeros(length(stats),1);MinorAxis=zeros(length(stats),1);
  for q = 1:length(stats);
    mi = m1==q;
    bb = stats(q).BoundingBox;
    i1 = round((bb(1)-1:bb(1)+bb(3)+1)); i1(i1<1) = []; i1(i1>size(chl555,2)) = []; 
    i2 = round((bb(2)-1:bb(2)+bb(4)+1)); i2(i2<1) = []; i2(i2>size(chl555,1)) = [];
    %imshow([m0(i2,i1) mi(i2,i1)]); drawnow; pause;
    M0{w} = m0(i2,i1); Mi{w} = mi(i2,i1);
    ABS0{w} = chl555(i2,i1); LABEL0{w} = chl488(i2,i1); 
    MajorAxis(w)=stats(q).MajorAxisLength;MinorAxis(w)=stats(q).MinorAxisLength;
    fr(w) = k;
    w = w + 1;
  end

end
  
save([fname,'_segmented.mat'],'ABS0','LABEL0','M0','Mi','fr','MajorAxis','MinorAxis')