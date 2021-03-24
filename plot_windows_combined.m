function plot_windows_combined(downsampled,dtime,dscore,depochs,sT,min_win,allPeakIds,conv_zt,path)

newFolder= fullfile(path,'proc_bins');

if exist(newFolder,'dir')
    path = newFolder;
else
    mkdir(newFolder);
    path = newFolder;
end



%get indexes for bins by default EEG  is 400 Hz
% fd = round(1/median(diff(dtime)));
%
% dffbinlen = 1:min_win*60*fd:numel(downsampled);
% eegbinlen =1:400*60*min_win:numel(sT.Time);
fs = 400; % for the EEG data
[b a] = butter(2,0.05);
downsampled = filtfilt(b,a,downsampled);
%
% peaksFilter = designfilt('lowpassiir','HalfPowerFrequency',0.2,'SampleRate',fd,'DesignMethod','butter','FilterOrder',12);
% downsampled = filtfilt(peaksFilter,downsampled);

%calculate it based on epochs, 60*min_win is the total number of seconds,
%each epoch is 5 seconds
step = (60*min_win)/5;
n_tenbin = fix(max(depochs)/step);
ep_idx = 0:step:n_tenbin*step;
for i=1:length(ep_idx)-1
    gidx(i,2) = find(depochs == ep_idx(i+1),1,'last');
    gidx(i,1) = find(depochs == ep_idx(i)+1,1);
    eidx(i,2) = find(sT.Epochs == ep_idx(i+1),1,'last');
    eidx(i,1) = find(sT.Epochs == ep_idx(i)+1,1);
end

dumb_arr = zeros(numel(downsampled,1));
dumb_arr(allPeakIds) = 1;

for i=1:size(gidx,1)
    clf
    tmp_g = downsampled(gidx(i,1):gidx(i,2));
    colors = dscore(gidx(i,1):gidx(i,2));
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
    sec_idx =find(diff(floor(dtime(gidx(i,1):gidx(i,2))))==1);
    xticks([0 ; sec_idx(100:100:500) ; length(tmp_g)]);
    xticklabels(num2cell([0:100:600]));
    xlabel('Time (sec)')
    grid off;
    xlim([0 length(tmp_g)]);
    set(gca,'FontName', 'Arial')
    hold all
    tmp_mark = dumb_arr(gidx(i,1):gidx(i,2));
    mark_idx = find(tmp_mark);
    plot(mark_idx,tmp_g(tmp_mark==1),'o','markerfacecolor',[0 0 0],'MarkerEdgeColor','None')
    pbaspect([1 0.25 0.25]);
    
    exportgraphics(gcf,fullfile(path,['Fiber-Window-bin-' num2str(i) '-ZT-' num2str(conv_zt(1)) '-to-' num2str(conv_zt(2)) '.pdf']),'Resolution',100,'ContentType','vector')
    clf
    temg = downsample(sT.EMG(eidx(i,1):eidx(i,2)),10);
    tscores = downsample(sT.Scores(eidx(i,1):eidx(i,2)),10);
    sec_idx =find(diff(downsample(sT.Time(eidx(i,1):eidx(i,2)),10))=='00:00:01');
    x = 1:length(temg);
    y = temg';
    z = zeros(size(x));
    lineColor = tscores';
    surface([x;x], [y;y], [z;z], [lineColor;lineColor],...
        'FaceColor', 'no',...
        'EdgeColor', 'interp',...
        'LineWidth', 1.5);
    grid off;
    caxis([1 4])
    xlim([0 length(temg)]);
    ylim([-500 500]);
    xticks([0 ; sec_idx(100:100:500) ;length(temg)]);
    xticklabels(num2cell([0:100:600]));
    xlabel('Time (sec)')
    set(gca,'FontName', 'Arial')
    ylabel('EMG (mV)');
    
    pbaspect([1 0.25 0.25]);
    exportgraphics(gcf,fullfile(path,['EMG-Window-bin-' num2str(i) '-ZT-' num2str(conv_zt(1)) '-to-' num2str(conv_zt(2)) '.pdf']),'Resolution',100,'ContentType','vector')
end

