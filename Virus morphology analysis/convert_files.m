function [ABS,LABEL] = convert_files(savename);

fnames = strcat(savename,'.tif');

for k = 1:5;
  ABS(:,:,k) = imread(fnames,2*(k-1)+1);
  LABEL(:,:,k) = imread(fnames,2*(k-1)+2);
end
savename=strcat(savename,'.mat');
save(savename,'ABS','LABEL');
