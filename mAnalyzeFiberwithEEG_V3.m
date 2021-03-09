clear
[gfile, path] = uigetfile('*.csv','Pick .csv file for GCAMP signal');
gdirpath = fullfile(path,gfile);
cd(path)
[efile, path] = uigetfile('*.txt', 'Pick your EEG file');
edirpath = fullfile(path,efile);
%  [efile, path] = uigetfile('*.edf', 'Pick your EEG.edf file');
%  edirpath = fullfile(path,efile);
[sfile, path] = uigetfile('*.txt', 'Pick your score file');
sdirpath = fullfile(path,sfile);


%% This is the ZT 0-3
ZT_bin = 1;
if ~ZT_bin
    T = readtable(edirpath);
end
%[hdr, T] = edfread(edirpath);   %dropping the edf for now, I don't think I
%can get the times to work correctly with the scores.
eeg_Date = T{:,1};
eeg_Time = T{:,2};
eeg_TimeStamp = T{:,3};
if ZT_bin
    [ZT0_ind] = find(eeg_Time=='21:00:00');
    [ZT3_ind] = find(eeg_Time=='00:00:00');
    disp('Processing ZT12 TO 15....')
else
    [ZT0_ind] = find(eeg_Time=='09:00:00');
    [ZT3_ind] = find(eeg_Time=='12:00:00');
    disp('Processing ZT0 TO 3...')
end
eeg_TTL = T{:,7};
[TTL_ind] = find(contains(eeg_TTL,'TTL: Rise;'));
start_ind = find(TTL_ind>ZT0_ind(1));
stop_ind = find(TTL_ind<ZT3_ind(end));
fiber_cut = start_ind(1)-1:stop_ind(end)+1; %number of TLL pulses to be extracted
subTTL_ind = TTL_ind(start_ind(1)-1:stop_ind(end)+1);
EEG_ZT0_3 = T(subTTL_ind(1):subTTL_ind(end),:);
eeg_sTimeStamp = EEG_ZT0_3.TimeStamp;

%load the gcamp data
if ~ZT_bin
    G = readtable(gdirpath);
end
gTTL = round(G{:,4});
gTTL_ind = find(diff(gTTL)>0); %this gives us the transitions from 0 to 1 except the first pulse
gTTL_ind = [1; gTTL_ind]; %we are adding 1 to the first pulse to correctly index the cutoffs from EEG TTL channel
fiber_ZT0_3 = G(gTTL_ind(fiber_cut(1))+1:gTTL_ind(fiber_cut(end))-1,:);

%Analyze the fiber signal to calculate Z score
raw_signal = fiber_ZT0_3.AnalogIn__Ch_1AIn_1_Dem_AOut_2__LowPass;
raw_reference = fiber_ZT0_3.AnalogIn__Ch_1AIn_1_Dem_AOut_1__LowPass;
reg = polyfit(raw_reference, raw_signal, 1);
a = reg(1);
b = reg(2);
controlFit = a.*raw_reference +b;
normDat = (raw_signal - controlFit)./ controlFit; %this gives deltaF/F
normDat = normDat * 100; % get %nce + b;

%load the scoring data
if ~ZT_bin
    S = readtable(sdirpath);
end
s_scores = S{:,5};
s_TimeStamp = S{:,3};


[inter_Stamp iEEG iscore] = intersect(eeg_sTimeStamp,s_TimeStamp);
combined_TimeScores = zeros(size(eeg_sTimeStamp,1),3);
combined_TimeScores(:,1) = eeg_sTimeStamp;

for i=1:length(iscore)-1
    epoch_start = find(eeg_sTimeStamp == s_TimeStamp(iscore(i)));
    epoch_end = find(eeg_sTimeStamp == s_TimeStamp(iscore(i)+1));
    scoring = ones(size(epoch_start(1):epoch_end(1)-1,2),1);
    combined_TimeScores(epoch_start(1):epoch_end(1)-1,2) = scoring*s_scores(iscore(i));
    combined_TimeScores(epoch_start(1):epoch_end(1)-1,3) = i;
end
ind =find(combined_TimeScores(:,2)>0);
newTable = EEG_ZT0_3(ind(1):ind(end),:);
newTable.Scores = combined_TimeScores(ind(1):ind(end),2);
newTable.Epochs = combined_TimeScores(ind(1):ind(end),3);


%newTable variable is the EEG data with the scores and epochs
% for the same duration we have length fiber data and length EEG data. a
% simple conversion to find the length of 1000 values in
% the EEG data would be multiplying this number by the lenght of the fiber
% data and dividing it by the lenght of EEG data.
clear cut_normDat
fiberstart = floor((ind(1)*length(normDat))/length(eeg_sTimeStamp));
fiberend   = floor((ind(end)*length(normDat))/length(eeg_sTimeStamp));
cut_normDat = normDat(fiberstart:fiberend);
%now find all the epochs
fiber_ind = floor((find(diff(newTable.Epochs)==1)*length(cut_normDat))/length(ind));
cut_normDat = [cut_normDat zeros(size(cut_normDat,1),1) zeros(size(cut_normDat,1),1)]; %use repmat
cut_normDat(1:fiber_ind(1),2) =1 ;
cut_normDat(1:fiber_ind(1),3) = s_scores(iscore(1));
fiber_ind = [fiber_ind ;length(cut_normDat)];
for i=1:length(fiber_ind)-1
    cut_normDat(fiber_ind(i)+1:fiber_ind(i+1),2) = i+1;
    cut_normDat(fiber_ind(i)+1:fiber_ind(i+1),3) = s_scores(iscore(i+1));
end

modFiber = fiber_ZT0_3(fiberstart:fiberend,:);
modFiber.normData = cut_normDat(:,1);
modFiber.Epochs = cut_normDat(:,2);
modFiber.Scores = cut_normDat(:,3);

ref_raw = modFiber.AnalogIn__Ch_1AIn_1_Dem_AOut_1__LowPass;
g_raw = modFiber.AnalogIn__Ch_1AIn_1_Dem_AOut_2__LowPass;
time = modFiber.Time_s_;
fs = 1/median(diff(time)); %calculates the sampling frequency

bleachingFilter = designfilt('lowpassiir', 'HalfPowerFrequency', 0.1, 'SampleRate', fs, 'DesignMethod', 'butter', 'FilterOrder', 12);
rLowpass = filtfilt(bleachingFilter, ref_raw);
sLowpass = filtfilt(bleachingFilter, g_raw);
rFit = fit(time, rLowpass, fittype('exp1'));
sFit = fit(time, sLowpass, fittype('exp1'));
rBleaching = rFit(time);
sBleaching = sFit(time);
rCorrected = ref_raw - rBleaching;
sCorrected = g_raw - sBleaching;

% Correct for movement artifacts.
f = sCorrected - rCorrected;

n0 = nanmin(round(600 * fs), numel(f));
f0 = movmean(f,n0);
n1 = nanmin(round(600 * fs), numel(f));
f1 = movstd(f,n1);
dff = (f - f0) ./ f1;
peaksFilter = designfilt('lowpassiir', 'HalfPowerFrequency', 0.2, 'SampleRate', fs, 'DesignMethod', 'butter', 'FilterOrder', 12);
peaksLowpass = filtfilt(peaksFilter, dff);
lowpassFilter = designfilt('lowpassiir', 'HalfPowerFrequency', 2, 'SampleRate', fs, 'DesignMethod', 'butter', 'FilterOrder', 12);
dffLowpass = filtfilt(lowpassFilter, dff);
modFiber.dff = dff;
modFiber.dffLowpass = dffLowpass;
modFiber.peaks = peaksLowpass;

bandpassFilter = designfilt('bandpassfir', 'CutoffFrequency1', 0.02, 'CutoffFrequency2', 0.2, 'SampleRate', fs, 'DesignMethod', 'window', 'FilterOrder', 12);
dffBandpass = filtfilt(bandpassFilter, dff);
modFiber.dffbandpass = dffBandpass;



%save processed files so no need to open and process from the original

if ZT_bin
    writetable(newTable, fullfile(path,'EEG_proc_ZT12_15_Sliced.txt'));
    writetable(modFiber, fullfile(path,'GCaMP_proc_ZT12_15_Sliced.txt'));
else
    writetable(newTable, fullfile(path,'EEG_proc_ZT0_3_Sliced.txt'));
    writetable(modFiber, fullfile(path,'GCaMP_proc_ZT0_3_Sliced.txt'));
end

%calculate transients eban-rothschild sytle
%each epoch is 5 seconds to reach 10 mins of data we need 120 epochs in
%total to reach 10 min of data
%find indices for 10 mins
clear gidx eidx
n_tenbin = fix(max(modFiber.Epochs)/120); %this around 18 which is 3 hours of 10 min bins.
ep_idx = 0:120:n_tenbin*120;
for i =1:length(ep_idx)-1
    gidx(i,2) =  find(modFiber.Epochs==ep_idx(i+1),1,'last');
    gidx(i,1) = find(modFiber.Epochs==ep_idx(i)+1,1);
    eidx(i,2) = find(newTable.Epochs==ep_idx(i+1),1,'last');
    eidx(i,1) = find(newTable.Epochs==ep_idx(i)+1,1);
end

[b a] = butter(2,0.05);

for i=1:size(gidx,1)
    clf
    subplot(311)
    tmp_g = downsample(modFiber.dffLowpass(gidx(i,1):gidx(i,2)),28); %brings it down to 40 Hz
    tmp_g = filtfilt(b,a,tmp_g);
    colors = downsample(modFiber.Scores(gidx(i,1):gidx(i,2)),28);
    cmap = brewermap(4,'Set1');
    colormap(cmap)
    x = 1:length(tmp_g);
    y = tmp_g';
    z = zeros(size(x));
    lineColor = colors';
    surface([x;x], [y;y], [z;z], [lineColor;lineColor],...
        'FaceColor', 'no',...
        'EdgeColor', 'interp',...
        'LineWidth', 1.5);
    caxis([1 4])
    sec_idx =find(diff(floor(downsample(modFiber.Time_s_(gidx(i,1):gidx(i,2)),28)))==1);
    xticks([0 ; sec_idx(100:100:500) ; length(tmp_g)]);
    xticklabels(num2cell([0:100:600]));
    xlabel('Time (sec)')
    grid off;
    xlim([0 length(tmp_g)]);
    set(gca,'FontName', 'Arial')
    subplot(312)
    teeg = newTable.EEG1(eidx(i,1):eidx(i,2));
    tscores = newTable.Scores(eidx(i,1):eidx(i,2));
    x = 1:length(teeg);
    y = teeg';
    z = zeros(size(x));
    lineColor = tscores';
    surface([x;x], [y;y], [z;z], [lineColor;lineColor],...
        'FaceColor', 'no',...
        'EdgeColor', 'interp',...
        'LineWidth', 1.5);
    grid off;
    caxis([1 4])
    xlim([0 length(teeg)])
    sec_idx =find(diff(newTable.Time(eidx(i,1):eidx(i,2)))=='00:00:01');
    xticks([0 ; sec_idx(100:100:500) ;length(teeg)]);
     xticklabels(num2cell([0:100:600]));
    xlabel('Time (sec)')
    set(gca,'FontName', 'Arial')
    subplot(313)
    teeg = newTable.EMG(eidx(i,1):eidx(i,2));
    tscores = newTable.Scores(eidx(i,1):eidx(i,2));
    x = 1:length(teeg);
    y = teeg';
    z = zeros(size(x));
    lineColor = tscores';
    surface([x;x], [y;y], [z;z], [lineColor;lineColor],...
        'FaceColor', 'no',...
        'EdgeColor', 'interp',...
        'LineWidth', 1.5);
    grid off;
    caxis([1 4])
    xlim([0 length(teeg)])
     xticks([0 ; sec_idx(100:100:500) ;length(teeg)]);
     xticklabels(num2cell([0:100:600]));
    xlabel('Time (sec)')
    set(gca,'FontName', 'Arial')
    if ~ZT_bin
    exportgraphics(gcf,fullfile(path,['zt0_3_10minWindows_bin_' num2str(i) '_' '.pdf']),'Resolution',100,'ContentType','vector')
    exportgraphics(gcf,fullfile(path,['zt0_3_10minWindows_bin_' num2str(i) '_' '.png']),'Resolution',300,'ContentType','image')
    else
        exportgraphics(gcf,fullfile(path,['zt12_15_10minWindows_bin_' num2str(i) '_' '.png']),'Resolution',100,'ContentType','image')
        exportgraphics(gcf,fullfile(path,['zt12_15_10minWindows_bin_' num2str(i) '_' '.pdf']),'Resolution',100,'ContentType','vector')
     end
end



 calc_transients = downsample(modFiber.dffLowpass,28); %brings it down to 40 Hz
 calc_transients = filtfilt(b,a,calc_transients);
 [peaks locs] = findpeaks(calc_transients,'MinPeakProminence',0.5,'MinPeakWidth',30);
 transient_scores = downsample(modFiber.Scores,28);
 transient_epochs = downsample(modFiber.Epochs,28);

 %find frequencies, calculate how many seconds in 3 hours for each epoch
 for i=1:4
     transients(i) = length(find(transient_scores(locs) == i));
     secs(i) = length(unique(transient_epochs(transient_scores==i)))*5;
     freqs(i) = transients(i)/secs(i);
 end
 
 if ~ZT_bin
     save('zt03_transient.mat','transients','secs','freqs')
 else
     save('zt12_15_transient.mat','transients','secs','freqs')
 end

%%
%plot 20 min windows for each



%
% subplot(311)
% g = scatter(1:length(downsample(modFiber.dff,10)),downsample(modFiber.dff,10),[],downsample(modFiber.Scores,10),'filled');
% xlim([0,length(downsample(modFiber.dff,10))]);
% g.SizeData = 2;
% downstime = downsample(modFiber.Time_s_,10);
% downstime = downstime - downstime(1);
% last =floor(downstime(end)/60);
% last = find(downstime<last*60,1,'last');
% interval = find(downstime<60*10,1,'last');
% tickArray = 0:interval:last;
% xticks(tickArray);
% tickLabels = 0:10:(length(tickArray)-1)*10;
% tickLabels = compose('%d',tickLabels);
% xticklabels(tickLabels);
% ylabel('df/f')
% xlabel('Time (min)')
% title(gfile)
% %legend('quiet wake','nonrem','rem','active wake')
% subplot(312)
% e = scatter(1:length(newTable.EEG1), newTable.EEG1,[],newTable.Scores,'filled');
% e.SizeData = 2;
% xlim([0, length(newTable.EEG1)]);
% time = newTable.TimeFromStart;
% time = time - time(1);
% last = (time(end)/60);
% last = find(time<last*60,1,'last');
% interval = find(time<60*10,1,'last');
% tickArray = 0:interval:last;
% xticks(tickArray);
% tickLabels = 0:10:(length(tickArray)-1)*10;
% tickLabels = compose('%d',tickLabels);
% xticklabels(tickLabels);
% ylabel('EEG 10 mV')
% xlabel('Time (min)')
% subplot(313)
% m = scatter(1:length(newTable.EMG), newTable.EMG,[],newTable.Scores,'filled');
% m.SizeData = 2;
% xlim([0, length(newTable.EMG)]);
% xticks(tickArray);
% xticklabels(tickLabels);
% ylabel('EMG 10 mV')
% xlabel('Time (min)')



%downsample it by a factor of 50, brings it down to 1250/50 to 25Hz
downsampled = downsample(modFiber.dffLowpass,50);
depochs = downsample(modFiber.Epochs,50);
dscore = downsample(modFiber.Scores,50);
%NonRem                 	2
%REM                    	3
%Unscored               	255
%Wake                   	1
%find where nrem is
[row,col] = find(dscore==2);
[drow ,dcol] = find(diff(row)>1);
%drow are the locations of NREM switches.
%check what comes before it to seperate
clear comb_ind
comb_ind.rem2nrem = [];
comb_ind.nrem2active = [];
comb_ind.nrem2rem = [];
comb_ind.active2nrem = [];
for i=1:length(drow)
    %after nrem
    if dscore(row(drow(i))-5) == 2 %leaving NREM
        if dscore(row(drow(i))+5) == 4 || dscore(row(drow(i))+5) == 1 %NREM to active wake
            comb_ind.nrem2active(end+1) = row(drow(i));
        elseif dscore(row(drow(i))+5) ==3 %NREM to REM
            comb_ind.nrem2rem(end+1) = row(drow(i));
        end
    end
    %before nrem
    if dscore(row(drow(i)+1)-5)== 3 %REM to NREM
        if dscore(row(drow(i)+1)+5) == 2 %internal check
            comb_ind.rem2nrem(end+1) = row(drow(i)+1);
        end
    elseif dscore(row(drow(i)+1)-5) ==4 || dscore(row(drow(i)+1)-5) ==1%active to NREM
        if dscore(row(drow(i)+1)+5) == 2 %internal check
            comb_ind.active2nrem(end+1) = row(drow(i)+1);
        end
    end
end


%add +/- 30 seconds
fn = fieldnames(comb_ind);
clear comb_dat comb_scores comb_depochs
for k=1:numel(fn)
    comb_dat.(fn{k}) = [];
    comb_scores.(fn{k}) = [];
    comb_depochs.(fn{k}) = [];
    %comb_dat(k).scores = [];
    %comb_dat(k).ident = fn{k};
    for i=1:length(comb_ind.(fn{k}))
        if comb_ind.(fn{k})(i)-600<0 || ~ comb_ind.(fn{k})(i)+600 > length(downsampled)
            continue
        else
            comb_dat.(fn{k})(end+1,:) = downsampled(comb_ind.(fn{k})(i)-600:comb_ind.(fn{k})(i)+600);
            comb_scores.(fn{k})(end+1,:) = dscore(comb_ind.(fn{k})(i)-600:comb_ind.(fn{k})(i)+600);
            comb_depochs.(fn{k})(end+1,:) = depochs(comb_ind.(fn{k})(i)-600:comb_ind.(fn{k})(i)+600);
        end
    end
end

if ZT_bin
    save([gdirpath(1:end-4) '_proc_transition_just_wake_ZT12_15.mat'],'comb_dat','comb_scores','comb_depochs')
    disp('saved ZT12-15')
else
    save([gdirpath(1:end-4) '_proc_transition_just_wake_ZT0_3.mat'],'comb_dat','comb_scores','comb_depochs')
    disp('saved ZT0-3')
end
%%







%Calculates the mean dff (and other values, column 9=dff in the table)
%based on scored state.
state_avgs = varfun(@mean,modFiber,'GroupingVariables','Scores');

%cutting modFiber into 1 hr bins
modfiber_hour = round(height(modFiber)/3);
modFiber0_1 = modFiber(1:modfiber_hour,:);
modFiber1_2 = modFiber(modfiber_hour:(modfiber_hour*2),:);
modFiber2_3 = modFiber((modfiber_hour*2):end,:);

%cutting eeg/emg (newTable) into 1 hr bins
EEG_hour = round(height(newTable)/3);
EEG_hr0_1 = newTable(1:EEG_hour,:);
EEG_hr1_2 = newTable(EEG_hour:(EEG_hour*2),:);
EEG_hr2_3 = newTable((EEG_hour*2):end,:);

%taking means of df/f by state
state_avgs0_1 = varfun(@mean,modFiber0_1,'GroupingVariables','Scores');
state_avgs1_2 = varfun(@mean,modFiber1_2,'GroupingVariables','Scores');
state_avgs2_3 = varfun(@mean,modFiber2_3,'GroupingVariables','Scores');
%collect these into one table
hourly_mean_dffs = outerjoin(state_avgs0_1(:,[1,9]), (outerjoin(state_avgs1_2(:,[1,9]), state_avgs2_3(:,[1,9]), 'Keys','Scores','MergeKeys',true)), 'Keys', 'Scores','MergeKeys',true);
hourly_mean_dffs.Properties.VariableNames{'mean_dff'} = 'mean_dff 0_1';
hourly_mean_dffs.Properties.VariableNames{'mean_dff_left'} = 'mean_dff 1_2';
hourly_mean_dffs.Properties.VariableNames{'mean_dff_right'} = 'mean_dff 2_3';

%taking medians of df/f by state
state_med_avgs0_1 = varfun(@median,modFiber0_1,'GroupingVariables','Scores');
state_med_avgs1_2 = varfun(@median,modFiber1_2,'GroupingVariables','Scores');
state_med_avgs2_3 = varfun(@median,modFiber2_3,'GroupingVariables','Scores');
%collect these into one table
hourly_med_dffs = outerjoin(state_med_avgs0_1(:,[1,9]), (outerjoin(state_med_avgs1_2(:,[1,9]), state_med_avgs2_3(:,[1,9]), 'Keys','Scores','MergeKeys',true)), 'Keys', 'Scores','MergeKeys',true);
hourly_med_dffs.Properties.VariableNames{'median_dff'} = 'median_dff 0_1';
hourly_med_dffs.Properties.VariableNames{'median_dff_left'} = 'median_dff 1_2';
hourly_med_dffs.Properties.VariableNames{'median_dff_right'} = 'median_dff 2_3';

%save as excel files
writetable(hourly_mean_dffs, strcat('/Users/benjaminbell/Desktop/Fiber Photometry/', strcat(sfile,' ',' 0 - 3 hourly means.xlsx')));
writetable(hourly_med_dffs, strcat('/Users/benjaminbell/Desktop/Fiber Photometry/', strcat(sfile,' ',' 0 - 3 hourly medians.xlsx')));

%plot each hour separately
subplot(311)
g = scatter(1:length(downsample(modFiber0_1.dff,10)),downsample(modFiber0_1.dff,10),[],downsample(modFiber0_1.Scores,10),'filled');
xlim([0,length(downsample(modFiber0_1.dff,10))]);
g.SizeData = 2;
downstime = downsample(modFiber0_1.Time_s_,10);
downstime = downstime - downstime(1);
last =floor(downstime(end)/60);
last = find(downstime<last*60,1,'last');
interval = find(downstime<60*10,1,'last');
tickArray = 0:interval:last;
xticks(tickArray);
tickLabels = 0:10:(length(tickArray)-1)*10;
tickLabels = compose('%d',tickLabels);
xticklabels(tickLabels);
ylabel('df/f')
xlabel('Time (min)')
title('0 - 1')

subplot(312)
g = scatter(1:length(downsample(modFiber1_2.dff,10)),downsample(modFiber1_2.dff,10),[],downsample(modFiber1_2.Scores,10),'filled');
xlim([0,length(downsample(modFiber1_2.dff,10))]);
g.SizeData = 2;
downstime = downsample(modFiber1_2.Time_s_,10);
downstime = downstime - downstime(1);
last =floor(downstime(end)/60);
last = find(downstime<last*60,1,'last');
interval = find(downstime<60*10,1,'last');
tickArray = 0:interval:last;
xticks(tickArray);
tickLabels = 0:10:(length(tickArray)-1)*10;
tickLabels = compose('%d',tickLabels);
xticklabels(tickLabels);
ylabel('df/f')
xlabel('Time (min)')
title('1 - 2')

subplot(313)
g = scatter(1:length(downsample(modFiber2_3.dff,10)),downsample(modFiber2_3.dff,10),[],downsample(modFiber2_3.Scores,10),'filled');
xlim([0,length(downsample(modFiber2_3.dff,10))]);
g.SizeData = 2;
downstime = downsample(modFiber2_3.Time_s_,10);
downstime = downstime - downstime(1);
last =floor(downstime(end)/60);
last = find(downstime<last*60,1,'last');
interval = find(downstime<60*10,1,'last');
tickArray = 0:interval:last;
xticks(tickArray);
tickLabels = 0:10:(length(tickArray)-1)*10;
tickLabels = compose('%d',tickLabels);
xticklabels(tickLabels);
ylabel('df/f')
xlabel('Time (min)')
title('2 - 3')
cb = colorbar;
cb.Limits = [1,4];
cb.Ticks = [1,2,3,4];
cb.Location = 'east';
onehrplots = gcf;

%save the hour length figures as .png
onehrplots = gcf;
exportgraphics(onehrplots,strcat('/Users/benjaminbell/Desktop/Fiber Photometry/', strcat(sfile,' ',' 0 - 3 dff traces.png')));
clf


%ploting 1 hr window of df/f, eeg, emg
%plotting 1hr of df/f, EEG, and EMG
subplot(311)
g = scatter(1:length(downsample(modFiber2_3.dff,10)),downsample(modFiber2_3.dff,10),[],downsample(modFiber2_3.Scores,10),'filled');
xlim([0,length(downsample(modFiber2_3.dff,10))]);
ylim([-2, 5]);
g.SizeData = 2;
downstime = downsample(modFiber2_3.Time_s_,10);
downstime = downstime - downstime(1);
last =floor(downstime(end)/60);
last = find(downstime<last*60,1,'last');
interval = find(downstime<60*10,1,'last');
tickArray = 0:interval:last;
xticks(tickArray);
tickLabels = 0:10:(length(tickArray)-1)*10;
tickLabels = compose('%d',tickLabels);
xticklabels(tickLabels);
ylabel('df/f')
xlabel('Time (min)')
title('2 - 3')
cb = colorbar;
cb.Limits = [1,4];
cb.Ticks = [1,2,3,4];
cb.Location = 'manual';
cb.Position = [.1 .1 0 0];


subplot(312)
e = scatter(1:length(EEG_hr2_3.EEG1), EEG_hr2_3.EEG1,[],EEG_hr2_3.Scores,'filled');
e.SizeData = 2;
xlim([0, length(EEG_hr2_3.EEG1)]);
time = EEG_hr2_3.TimeFromStart;
time = time - time(1);
last = (time(end)/60);
last = find(time<last*60,1,'last');
interval = find(time<60*10,1,'last');
tickArray = 0:interval:last;
xticks(tickArray);
tickLabels = 0:10:(length(tickArray)-1)*10;
tickLabels = compose('%d',tickLabels);
xticklabels(tickLabels);
ylabel('EEG 10 mV')
xlabel('Time (min)')
subplot(313)
m = scatter(1:length(EEG_hr2_3.EMG), EEG_hr2_3.EMG,[],EEG_hr2_3.Scores,'filled');
m.SizeData = 2;
xlim([0, length(EEG_hr2_3.EMG)]);
xticks(tickArray);
xticklabels(tickLabels);
ylabel('EMG 10 mV')
xlabel('Time (min)')

clf
%% Now ZT 12-15

T = readtable(edirpath);
%[hdr, T] = edfread(edirpath);   %dropping the edf for now, I don't think I
%can get the times to work correctly with the scores.
eeg_Date = T{:,1};
eeg_Time = T{:,2};
eeg_TimeStamp = T{:,3};
[ZT0_ind] = find(eeg_Time=='21:00:00');
[ZT3_ind] = find(eeg_Time=='00:00:00');
eeg_TTL = T{:,7};
[TTL_ind] = find(contains(eeg_TTL,'TTL: Rise;'));
start_ind = find(TTL_ind>ZT0_ind(1));
stop_ind = find(TTL_ind<ZT3_ind(end));
fiber_cut = start_ind(1)-1:stop_ind(end)+1; %number of TLL pulses to be extracted
subTTL_ind = TTL_ind(start_ind(1)-1:stop_ind(end)+1);
EEG_ZT0_3 = T(subTTL_ind(1):subTTL_ind(end),:);
eeg_sTimeStamp = EEG_ZT0_3.TimeStamp;

%load the gcamp data
G = readtable(gdirpath);
gTTL = round(G{:,4});
gTTL_ind = find(diff(gTTL)>0); %this gives us the transitions from 0 to 1 except the first pulse
gTTL_ind = [1; gTTL_ind]; %we are adding 1 to the first pulse to correctly index the cutoffs from EEG TTL channel
fiber_ZT0_3 = G(gTTL_ind(fiber_cut(1))+1:gTTL_ind(fiber_cut(end))-1,:);

%Analyze the fiber signal to calculate Z score
raw_signal = fiber_ZT0_3.AnalogIn__Ch_1AIn_1_Dem_AOut_2__LowPass;
raw_reference = fiber_ZT0_3.AnalogIn__Ch_1AIn_1_Dem_AOut_1__LowPass;
reg = polyfit(raw_reference, raw_signal, 1);
a = reg(1);
b = reg(2);
controlFit = a.*raw_reference +b;
normDat = (raw_signal - controlFit)./ controlFit; %this gives deltaF/F
normDat = normDat * 100; % get %nce + b;

%load the scoring data
S = readtable(sdirpath);
s_scores = S{:,5};
s_TimeStamp = S{:,3};


[inter_Stamp iEEG iscore] = intersect(eeg_sTimeStamp,s_TimeStamp);
combined_TimeScores = zeros(size(eeg_sTimeStamp,1),3);
combined_TimeScores(:,1) = eeg_sTimeStamp;

for i=1:length(iscore)-1
    epoch_start = find(eeg_sTimeStamp == s_TimeStamp(iscore(i)));
    epoch_end = find(eeg_sTimeStamp == s_TimeStamp(iscore(i)+1));
    scoring = ones(size(epoch_start(1):epoch_end(1)-1,2),1);
    combined_TimeScores(epoch_start(1):epoch_end(1)-1,2) = scoring*s_scores(iscore(i));
    combined_TimeScores(epoch_start(1):epoch_end(1)-1,3) = i;
end
ind =find(combined_TimeScores(:,2)>0);
newTable = EEG_ZT0_3(ind(1):ind(end),:);
newTable.Scores = combined_TimeScores(ind(1):ind(end),2);
newTable.Epochs = combined_TimeScores(ind(1):ind(end),3);


%newTable variable is the EEG data with the scores and epochs
% for the same duration we have length fiber data and length EEG data. a
% simple conversion to find the length of 1000 values in
% the EEG data would be multiplying this number by the lenght of the fiber
% data and dividing it by the lenght of EEG data.
clear cut_normDat
fiberstart = floor((ind(1)*length(normDat))/length(eeg_sTimeStamp));
fiberend   = floor((ind(end)*length(normDat))/length(eeg_sTimeStamp));
cut_normDat = normDat(fiberstart:fiberend);
%now find all the epochs
fiber_ind = floor((find(diff(newTable.Epochs)==1)*length(cut_normDat))/length(ind));
cut_normDat = [cut_normDat zeros(size(cut_normDat,1),1) zeros(size(cut_normDat,1),1)]; %use repmat
cut_normDat(1:fiber_ind(1),2) =1 ;
cut_normDat(1:fiber_ind(1),3) = s_scores(iscore(1));
fiber_ind = [fiber_ind ;length(cut_normDat)];
for i=1:length(fiber_ind)-1
    cut_normDat(fiber_ind(i)+1:fiber_ind(i+1),2) = i+1;
    cut_normDat(fiber_ind(i)+1:fiber_ind(i+1),3) = s_scores(iscore(i+1));
end

modFiber = fiber_ZT0_3(fiberstart:fiberend,:);
modFiber.normData = cut_normDat(:,1);
modFiber.Epochs = cut_normDat(:,2);
modFiber.Scores = cut_normDat(:,3);

ref_raw = modFiber.AnalogIn__Ch_1AIn_1_Dem_AOut_1__LowPass;
g_raw = modFiber.AnalogIn__Ch_1AIn_1_Dem_AOut_2__LowPass;
time = modFiber.Time_s_;
fs = 1/median(diff(time)); %calculates the sampling frequency

bleachingFilter = designfilt('lowpassiir', 'HalfPowerFrequency', 0.1, 'SampleRate', fs, 'DesignMethod', 'butter', 'FilterOrder', 12);
rLowpass = filtfilt(bleachingFilter, ref_raw);
sLowpass = filtfilt(bleachingFilter, g_raw);
rFit = fit(time, rLowpass, fittype('exp1'));
sFit = fit(time, sLowpass, fittype('exp1'));
rBleaching = rFit(time);
sBleaching = sFit(time);
rCorrected = ref_raw - rBleaching;
sCorrected = g_raw - sBleaching;

% Correct for movement artifacts.
f = sCorrected - rCorrected;

n0 = nanmin(round(600 * fs), numel(f));
f0 = movmean(f,n0);
n1 = nanmin(round(600 * fs), numel(f));
f1 = movstd(f,n1);
dff = (f - f0) ./ f1;
peaksFilter = designfilt('lowpassiir', 'HalfPowerFrequency', 0.2, 'SampleRate', fs, 'DesignMethod', 'butter', 'FilterOrder', 12);
peaksLowpass = filtfilt(peaksFilter, dff);
lowpassFilter = designfilt('lowpassiir', 'HalfPowerFrequency', 2, 'SampleRate', fs, 'DesignMethod', 'butter', 'FilterOrder', 12);
dffLowpass = filtfilt(lowpassFilter, dff);
modFiber.dff = dff;
modFiber.dffLowpass = dffLowpass;
modFiber.peaks = peaksLowpass;

bandpassFilter = designfilt('bandpassfir', 'CutoffFrequency1', 0.02, 'CutoffFrequency2', 0.2, 'SampleRate', fs, 'DesignMethod', 'window', 'FilterOrder', 12);
dffBandpass = filtfilt(bandpassFilter, dff);
modFiber.dffbandpass = dffBandpass;

%
% subplot(311)
% g = scatter(1:length(downsample(modFiber.dff,10)),downsample(modFiber.dff,10),[],downsample(modFiber.Scores,10),'filled');
% xlim([0,length(downsample(modFiber.dff,10))]);
% g.SizeData = 2;
% downstime = downsample(modFiber.Time_s_,10);
% downstime = downstime - downstime(1);
% last =floor(downstime(end)/60);
% last = find(downstime<last*60,1,'last');
% interval = find(downstime<60*10,1,'last');
% tickArray = 0:interval:last;
% xticks(tickArray);
% tickLabels = 0:10:(length(tickArray)-1)*10;
% tickLabels = compose('%d',tickLabels);
% xticklabels(tickLabels);
% ylabel('df/f')
% xlabel('Time (min)')
% title(gfile)
% %legend('quiet wake','nonrem','rem','active wake')
% subplot(312)
% e = scatter(1:length(newTable.EEG1), newTable.EEG1,[],newTable.Scores,'filled');
% e.SizeData = 2;
% xlim([0, length(newTable.EEG1)]);
% time = newTable.TimeFromStart;
% time = time - time(1);
% last = (time(end)/60);
% last = find(time<last*60,1,'last');
% interval = find(time<60*10,1,'last');
% tickArray = 0:interval:last;
% xticks(tickArray);
% tickLabels = 0:10:(length(tickArray)-1)*10;
% tickLabels = compose('%d',tickLabels);
% xticklabels(tickLabels);
% ylabel('EEG 10 mV')
% xlabel('Time (min)')
% subplot(313)
% m = scatter(1:length(newTable.EMG), newTable.EMG,[],newTable.Scores,'filled');
% m.SizeData = 2;
% xlim([0, length(newTable.EMG)]);
% xticks(tickArray);
% xticklabels(tickLabels);
% ylabel('EMG 10 mV')
% xlabel('Time (min)')

%Calculates the mean dff (and other values, column 9=dff in the table)
%based on scored state.
state_avgs = varfun(@mean,modFiber,'GroupingVariables','Scores');

%cutting modFiber into 1 hr bins
modfiber_hour = round(height(modFiber)/3);
modFiber12_13 = modFiber(1:modfiber_hour,:);
modFiber13_14 = modFiber(modfiber_hour:(modfiber_hour*2),:);
modFiber14_15 = modFiber((modfiber_hour*2):end,:);

%cutting eeg/emg (newTable) into 1 hr bins
EEG_hour = round(height(newTable)/3);
EEG_hr12_13 = newTable(1:EEG_hour,:);
EEG_hr13_14 = newTable(EEG_hour:(EEG_hour*2),:);
EEG_hr14_15 = newTable((EEG_hour*2):end,:);

%taking means of df/f by state
state_avgs12_13 = varfun(@mean,modFiber12_13,'GroupingVariables','Scores');
state_avgs13_14 = varfun(@mean,modFiber13_14,'GroupingVariables','Scores');
state_avgs14_15 = varfun(@mean,modFiber14_15,'GroupingVariables','Scores');
%collect these into one table
hourly_mean_dffs = outerjoin(state_avgs12_13(:,[1,9]), (outerjoin(state_avgs13_14(:,[1,9]), state_avgs14_15(:,[1,9]), 'Keys','Scores','MergeKeys',true)), 'Keys', 'Scores','MergeKeys',true);
hourly_mean_dffs.Properties.VariableNames{'mean_dff'} = 'mean_dff 12_13';
hourly_mean_dffs.Properties.VariableNames{'mean_dff_left'} = 'mean_dff 13_14';
hourly_mean_dffs.Properties.VariableNames{'mean_dff_right'} = 'mean_dff 14_15';

%taking medians of df/f by state
state_med_avgs12_13 = varfun(@median,modFiber12_13,'GroupingVariables','Scores');
state_med_avgs13_14 = varfun(@median,modFiber13_14,'GroupingVariables','Scores');
state_med_avgs14_15 = varfun(@median,modFiber14_15,'GroupingVariables','Scores');
%collect these into one table
hourly_med_dffs = outerjoin(state_med_avgs12_13(:,[1,9]), (outerjoin(state_med_avgs13_14(:,[1,9]), state_med_avgs14_15(:,[1,9]), 'Keys','Scores','MergeKeys',true)), 'Keys', 'Scores','MergeKeys',true);
hourly_med_dffs.Properties.VariableNames{'median_dff'} = 'median_dff 12_13';
hourly_med_dffs.Properties.VariableNames{'median_dff_left'} = 'median_dff 13_14';
hourly_med_dffs.Properties.VariableNames{'median_dff_right'} = 'median_dff 14_15';

%save as excel files
writetable(hourly_mean_dffs, strcat('/Users/benjaminbell/Desktop/Fiber Photometry/', strcat(sfile,' ',' 12 - 15 hourly means.xlsx')));
writetable(hourly_med_dffs, strcat('/Users/benjaminbell/Desktop/Fiber Photometry/', strcat(sfile,' ',' 12 - 15 hourly medians.xlsx')));

%plot each hour separately
subplot(311)
g = scatter(1:length(downsample(modFiber12_13.dff,10)),downsample(modFiber12_13.dff,10),[],downsample(modFiber12_13.Scores,10),'filled');
xlim([0,length(downsample(modFiber12_13.dff,10))]);
g.SizeData = 2;
downstime = downsample(modFiber12_13.Time_s_,10);
downstime = downstime - downstime(1);
last =floor(downstime(end)/60);
last = find(downstime<last*60,1,'last');
interval = find(downstime<60*10,1,'last');
tickArray = 0:interval:last;
xticks(tickArray);
tickLabels = 0:10:(length(tickArray)-1)*10;
tickLabels = compose('%d',tickLabels);
xticklabels(tickLabels);
ylabel('df/f')
xlabel('Time (min)')
title('12 - 13')

subplot(312)
g = scatter(1:length(downsample(modFiber13_14.dff,10)),downsample(modFiber13_14.dff,10),[],downsample(modFiber13_14.Scores,10),'filled');
xlim([0,length(downsample(modFiber13_14.dff,10))]);
g.SizeData = 2;
downstime = downsample(modFiber13_14.Time_s_,10);
downstime = downstime - downstime(1);
last =floor(downstime(end)/60);
last = find(downstime<last*60,1,'last');
interval = find(downstime<60*10,1,'last');
tickArray = 0:interval:last;
xticks(tickArray);
tickLabels = 0:10:(length(tickArray)-1)*10;
tickLabels = compose('%d',tickLabels);
xticklabels(tickLabels);
ylabel('df/f')
xlabel('Time (min)')
title('13 - 14')

subplot(313)
g = scatter(1:length(downsample(modFiber14_15.dff,10)),downsample(modFiber14_15.dff,10),[],downsample(modFiber14_15.Scores,10),'filled');
xlim([0,length(downsample(modFiber14_15.dff,10))]);
g.SizeData = 2;
downstime = downsample(modFiber14_15.Time_s_,10);
downstime = downstime - downstime(1);
last =floor(downstime(end)/60);
last = find(downstime<last*60,1,'last');
interval = find(downstime<60*10,1,'last');
tickArray = 0:interval:last;
xticks(tickArray);
tickLabels = 0:10:(length(tickArray)-1)*10;
tickLabels = compose('%d',tickLabels);
xticklabels(tickLabels);
ylabel('df/f')
xlabel('Time (min)')
title('14 - 15')
cb = colorbar;
cb.Limits = [1,4];
cb.Ticks = [1,2,3,4];
cb.Location = 'east';
onehrplots = gcf;

%save the hour length figures as .png
onehrplots = gcf;
exportgraphics(onehrplots,strcat('/Users/benjaminbell/Desktop/Fiber Photometry/', strcat(sfile,' ',' 12 - 15 dff traces.png')));
%clf


%plotting 1hr of df/f, EEG, and EMG
subplot(311)
g = scatter(1:length(downsample(modFiber14_15.dff,10)),downsample(modFiber14_15.dff,10),[],downsample(modFiber14_15.Scores,10),'filled');
xlim([0,length(downsample(modFiber14_15.dff,10))]);
ylim([-2, 5]);
g.SizeData = 2;
downstime = downsample(modFiber14_15.Time_s_,10);
downstime = downstime - downstime(1);
last =floor(downstime(end)/60);
last = find(downstime<last*60,1,'last');
interval = find(downstime<60*10,1,'last');
tickArray = 0:interval:last;
xticks(tickArray);
tickLabels = 0:10:(length(tickArray)-1)*10;
tickLabels = compose('%d',tickLabels);
xticklabels(tickLabels);
ylabel('df/f')
xlabel('Time (min)')
title('14 - 15')
cb = colorbar;
cb.Limits = [1,4];
cb.Ticks = [1,2,3,4];
cb.Location = 'manual';
cb.Position = [.1 .1 0 0];


subplot(312)
e = scatter(1:length(EEG_hr14_15.EEG1), EEG_hr14_15.EEG1,[],EEG_hr14_15.Scores,'filled');
e.SizeData = 2;
xlim([0, length(EEG_hr14_15.EEG1)]);
time = EEG_hr14_15.TimeFromStart;
time = time - time(1);
last = (time(end)/60);
last = find(time<last*60,1,'last');
interval = find(time<60*10,1,'last');
tickArray = 0:interval:last;
xticks(tickArray);
tickLabels = 0:10:(length(tickArray)-1)*10;
tickLabels = compose('%d',tickLabels);
xticklabels(tickLabels);
ylabel('EEG 10 mV')
xlabel('Time (min)')
subplot(313)
m = scatter(1:length(EEG_hr14_15.EMG), EEG_hr14_15.EMG,[],EEG_hr14_15.Scores,'filled');
m.SizeData = 2;
xlim([0, length(EEG_hr14_15.EMG)]);
xticks(tickArray);
xticklabels(tickLabels);
ylabel('EMG 10 mV')
xlabel('Time (min)')

onehr_all3 = gcf;
%save the hour length figures as .pdf   *** this isn't working, just use
%the export buttons in the figure
% exportgraphics(onehr_all3,strcat('/Users/benjaminbell/Desktop/Fiber Photometry/', '120520 zt14-15 traces.pdf'),'ContentType', 'vector');
exportgraphics(onehr_all3, 'm6_120520 zt14-15 traces.pdf','ContentType', 'vector');

%clf