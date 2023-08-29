%Put files into the TIFF folders where the files locate
function process_data()
fnames = ls('*.tif');
fnames_split=split(fnames);
fnames_char=char(fnames_split);
dim=size(fnames_char);
num_files=dim(1);
%disp(num_files);
%for n=1
for n=1:num_files-1
    s=fnames_char(n,:);
    savename=s(1:strfind(s,'.tif')-1);
    [ABS,LABEL] = convert_files(savename);
    segment_virions(savename);
    [f0,a0] = extract_data(savename);
end
%figure
%scatter(f0,a0,'.');
%set(gca,'xscale','log','yscale','log');drawnow
%filename=strcat('Rep',num2str(m),'_',conc,'nM.fig');
%savefig(filename);
end


