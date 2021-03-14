function [rsfG] = export_proc(sfG,sT,path,frequency,input_points)

fs = 1/median(diff(sfG.Time_s_)); %sampling frequency

[p ,q] = rat(frequency/fs); 
rsfG.dff = resample(sfG.dff,sfG.Time_s_,frequency,p,q);
rsfG.Epoch = resample(sfG.Epochs,sfG.Time_s_,frequency,p,q);
rsfG.Scores = resample(sfG.Scores,sfG.Time_s_,frequency,p,q);

%estimate the zt time
time_arr= 0:23;
zttime = [15:23 0:14];
for i=1:numel(input_points)
conv_zt(i) = zttime(time_arr == str2num(input_points{i}(1:2)));
end

 writetable(newTable, fullfile(path,'EEG_proc_ZT12_15_Sliced.txt'));



