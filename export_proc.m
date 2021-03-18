function [rsfG] = export_proc(rsfG,sT,path,frequency,input_points)



%estimate the zt time
time_arr= 0:23;
zttime = [15:23 0:14];
for i=1:numel(input_points)
    conv_zt(i) = zttime(time_arr == str2num(input_points{i}(1:2)));
end



%resampling can introduce artifacts, fix those
found = 1;
n=2;
while found
    found = isempty(find(round(diff(rsfG.time(1:n)))==0,1));
    n=n+1;
end
rsfG(1:n,:) = [];
sT(1:round(400*n/100),:) = []; 


%save the processed files
writetable(rsfG, fullfile(path,['GCaMP_dff' 'ZT-' num2str(conv_zt(1)) '-to-' num2str(conv_zt(2)) '-at-' num2str(frequency) 'Hz-' 'processed.csv']));
writetable(sT, fullfile(path,['EEG-EMG-' 'ZT-' num2str(conv_zt(1)) '-to-' num2str(conv_zt(2)) '-at-400'  'Hz-' 'processed.csv']));

end
