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