function [f0,a0] = extract_data(fname);

load([fname,'_segmented.mat']);
K = length(ABS0);

for k = 1:K
  f = ABS0{k}; a = LABEL0{k}; m0 = M0{k}; mi = Mi{k};
  temp = f.*(1-m0); temp = temp(temp>0); f_md = median(temp(:)); f_mn = mean(temp(:));
  temp = a.*(1-m0); temp = temp(temp>0); a_md = median(temp(:)); a_mn = mean(temp(:));
  %temp = c.*(1-m0); temp = temp(temp>0); c_md = median(temp(:)); c_mn = mean(temp(:));
  %X = [100*f 100*a 2^16*m0];
  %imshow(uint16(X)); drawnow; pause;
  f0(k) = sum(sum((f-f_md).*mi));
  a0(k) = sum(sum((a-a_md).*mi));
 %c0(k) = sum(sum((c-c_md).*mi));
end