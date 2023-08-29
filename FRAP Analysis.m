load('ROI.mat') %Load ROI for bleached spot
fnames = {dir('*.tif').name};
dim=size(fnames);
num_files=dim(2);
numframes = 15;
roi_sum = zeros(num_files,numframes);
medianXroi = zeros(num_files,numframes);
sum_chart = zeros(num_files,1);
HAroi_sum = zeros(num_files,numframes);
HAmedianXroi = zeros(num_files,numframes);
sum_chartHA = zeros(num_files,1);

for i=1:num_files
    for j=1:numframes
        curframe = imread(fnames{i},j);
        sumintensity = sum(sum(double(curframe).*ROI)); %sum fluorescent activity within ROI
        medianXroi(i,j) = median(double(curframe).*~ROI,"all");%median fluorescent activity outside ROI
        roi_sum(i,j) = sumintensity;
        if j==1
            sum_chart(i)=sumintensity;
        end
    end
end

medianXroi_norm = medianXroi./(ones(num_files,15).*medianXroi(:,1)); %find percentage of photobleaching in each frame in comparison to frame 1
roi_sum_norm = (roi_sum_norm-(ones(num_files,15).*roi_sum_norm(:,2)))./(ones(num_files,15).*(roi_sum_norm(:,1)-roi_sum_norm(:,2))); % find percentage recovery in each frame after frame 1


Recovery = roi_sum_norm(:,15);
HA= sum_chart(:,1)-ones(num_files,1)*min(sum_chart(:,1));
HA_norm=HA./max(HA(:,1));
figure(2);
scatter(HA_norm, roi_sum_norm(:,15));
xlim([0 1]);
ylim([0 1]);
xlabel('Normalized HA Expression Level')
ylabel('Percentage recovery')
title('Fluorescence Recovery after Photobleaching (CR9114IgG)')
