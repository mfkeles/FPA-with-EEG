%Analyze Fiber Photometry Data Combined with EEG
[gfile, path] = uigetfile('*.csv','Pick .csv file for GCAMP signal');
gdirpath = fullfile(path,gfile);
cd(path)
[efile, path] = uigetfile('*.txt', 'Pick your EEG file');
edirpath = fullfile(path,efile);
[sfile, path] = uigetfile('*.txt', 'Pick your score file');
sdirpath = fullfile(path,sfile);

%load all the data, this takes some time
G = readtable(gdirpath);
T = readtable(edirpath);
T.Properties.VariableNames{7} ='TTL';
opts = detectImportOptions(sdirpath,'NumHeaderLines',10);
S = readtable(sdirpath,opts,'ReadVariableNames',true);
S.ExtraVar1= [];
sprintf(['Gcamp file: ' gfile '\n' 'EEG file: ' efile '\n' 'Scoring file: ' sfile])

%input the timepoints to be processed, 1st one is start 2nd one is end:
time_points = {{'21:00:00','00:00:00'},{'09:00:00','12:00:00'}};
normalizationoption=4; %4 seems to be working best, (f - fmedian)/f(std). user can check other options
frequency = 100;  %target frequency for fiber photometry data, this is crucial as diff animals might have slightly different freq
window = 20; %+/- seconds to be sliced for transitions
ds = 4; %downsampling factor (4 brings it down to 25Hz)
min_win = 10; %number of minutes to plot as bins
for i=1:numel(time_points)
    input_points = time_points{i};
    [sT, rsG] = slice_timepoints(input_points,T,G,S,frequency); %Align all EEG+EMG data with Fiber
    [rsfG] = calc_dff(rsG,normalizationoption); %normalize the data
    [rsfG,conv_zt] = export_proc(rsfG,sT,path,frequency,input_points); %saves the processed file
    [M,downsampled,dscore,dtime,depochs] = calc_transitions(rsfG,ds,window,conv_zt,path); %find and save all transitions
    [allPeakIds, uniquePeakIds] = calc_transients(M,downsampled,dscore,dtime,conv_zt,path);%calculate transient rates
    plot_windows_combined(downsampled,dtime,dscore,depochs,sT,min_win,allPeakIds,conv_zt,path) %plots 10 minute windows 
end