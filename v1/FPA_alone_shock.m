before_shock  = 10; %how many seconds to be sliced before the shock
after_shock = 20; %how many seconds to be sliced after the shock


%load the file
[gfile, path] = uigetfile('*.csv','Pick .csv file for GCAMP signal');
gdirpath = fullfile(path,gfile);
G = readtable(gdirpath);

%clear out the NaN values
G = rmmissing(G);
shocks = G.AIn_2;

%calculate the frequency
fs = round( 1/median(diff(G.Time_s_)));

%remove decay
bleachingFilter = designfilt('lowpassiir', 'HalfPowerFrequency', 0.1, 'SampleRate', fs, 'DesignMethod', 'butter', 'FilterOrder', 12);
rLowpass = filtfilt(bleachingFilter, G.AIn_1_Dem_AOut_1_);
sLowpass = filtfilt(bleachingFilter, G.AIn_1_Dem_AOut_2_);
rFit = fit(G.Time_s_, rLowpass, fittype('exp2'));
sFit = fit(G.Time_s_, sLowpass, fittype('exp2'));
rBleaching = rFit(G.Time_s_);
sBleaching = sFit(G.Time_s_);
rCorrected = G.AIn_1_Dem_AOut_1_ - rBleaching;
sCorrected = G.AIn_1_Dem_AOut_2_ - sBleaching;

%fit reference to the signal
fitdata = fit(rCorrected,sCorrected,fittype('poly1'),'Robust','on');
rCorrected = fitdata(rCorrected);

% Correct for movement artifacts.
f = sCorrected - rCorrected;

%lowpassfilter
lowpassFilter = designfilt('lowpassiir', 'HalfPowerFrequency', 4, 'SampleRate', fs, 'DesignMethod', 'butter', 'FilterOrder', 12);
fSmooth = filtfilt(lowpassFilter, f);


%calculate z-score
f0 = median(fSmooth);
f1 = std(fSmooth);
dff = (fSmooth - f0)./f1;

%calculate where the foot shocks are happening
invertedM = max(movmean(shocks,400)) - movmean(shocks,400);
[pks, locs] = findpeaks(invertedM,'MinPeakHeight',1,'MinPeakDistance',800);
idx= [];
for i=1:numel(locs)
    start = find(round(diff(shocks(locs(i)-1500:locs(i)+1500)))~=0,1);
    stop  = find(round(diff(shocks(locs(i)-1500:locs(i)+1500)))~=0,1,'last');
    idx(i,:) = locs(i) + [start stop] -1500;
end
shock_len = round(median(idx(:,2)-idx(:,1)));
idx(:,2) = idx(:,1) + shock_len;

%slice the before and after the foot shock
dff_shocks = [];
trial_size = [0; zeros(before_shock*fs,1) ;ones(shock_len,1); zeros(after_shock*fs,1)];
for i=1:size(idx,1)
    dff_shocks(i,:) = dff(idx(i,1) - before_shock*fs:idx(i,2) + after_shock*fs); 
end


miny = floor(min(min(dff_shocks)));
maxy = ceil(max(max(dff_shocks)));
cmap = brewermap(size(idx,1),'Spectral');
for i=1:size(idx,1)
    plot(dff_shocks(i,:),'Color',cmap(i,:),'LineWidth',1);
    x=[find(trial_size==1,1),find(trial_size==1,1,'last'),find(trial_size==1,1,'last'),find(trial_size==1,1)];
    y=[miny,miny,maxy,maxy];
    patch(x,y,'black','FaceAlpha',.1,'EdgeColor','none')
    ylim([miny maxy]);
    ylabel('Z-Score')
    xlabel('Time (s)')
    title(['Shock Number: ' num2str(i)]);
    xticks([0:fs*2:max(xlim)]);
    labelarr = cellfun(@num2str,(num2cell(0:2:max(xlim)/fs)),'un',0);
    xticklabels(labelarr);
    exportgraphics(gcf,fullfile(path,[gfile(1:end-4) '_shock_' num2str(i) '.pdf']),'Resolution',300,'ContentType','vector');
end

clf

for i=1:size(idx,1)
    plot(dff_shocks(i,:),'Color',cmap(i,:),'LineWidth',1);
    hold all
end
x=[find(trial_size==1,1),find(trial_size==1,1,'last'),find(trial_size==1,1,'last'),find(trial_size==1,1)];
y=[miny,miny,maxy,maxy];
patch(x,y,'black','FaceAlpha',.1,'EdgeColor','none')
ylim([miny maxy]);
ylabel('Z-Score')
xlabel('Time (s)')
title(['Shocks red to blue']);
xticks([0:fs*2:max(xlim)]);
labelarr = cellfun(@num2str,(num2cell(0:2:max(xlim)/fs)),'un',0);
xticklabels(labelarr);
exportgraphics(gcf,fullfile(path,[gfile(1:end-4) '_all_shocks_' num2str(i) '.pdf']),'Resolution',300,'ContentType','vector');

max_f = max(dff_shocks(:,find(trial_size==1,1):find(trial_size==1,1,'last')+fs*4),[],2);
exp_tab = array2table(max_f);
exp_tab.Shock_Numbers = [1:6]';

writetable(exp_tab, fullfile(path,[gfile(1:end-4) '_shock_maximum_values.csv']));



