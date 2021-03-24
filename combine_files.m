%combine processed data across animals
dirpath =  uigetdir('Pick folder containing animals');

%current directory to the
cd(dirpath)

%find all subdirectories
animalFolders = find_folders(dirpath);
n=1;
%within each mouse find the correct subdirectories
for i=1:numel(animalFolders)
    m=1;
    %enter each subFolder and find session folders
    sessionFolders = find_folders(fullfile(animalFolders(i).folder,animalFolders(i).name));
    %for each session folder identify files to be loaded
    for j=1:numel(sessionFolders)
        file_list = dir(fullfile(sessionFolders(j).folder,sessionFolders(j).name,'GCaMP_dffZT*'));
        if ~isempty(file_list) %then this is already processed
            mouse(n).session{m} = fullfile(sessionFolders(j).folder,sessionFolders(j).name);
            m=m+1;
        end
    end
    n=n+1;
    m=m+1;
end

%% calculate max dff and area under the curve
check = 1;
for i=1:numel(mouse)
    for j=1:numel(mouse(i).session)
        mean_paths = dir(fullfile(mouse(i).session{j},'meanDFF*'));
        area_paths = dir(fullfile(mouse(i).session{j},'area-ZT*'));
        %load all the data
        for ii=1:numel(mean_paths)
            %get ZTs
            B1 = regexp(mean_paths(ii).name,'\d*','Match');
            B2 = regexp(area_paths(ii).name,'\d*','Match');
            ind_mean = readtable(fullfile(mean_paths(ii).folder,mean_paths(ii).name));
            ind_area = readtable(fullfile(area_paths(ii).folder,area_paths(ii).name));
            if check
                ind_mean.mouse = i;
                ind_mean.session= j;
                ind_mean.ZT = {['ZT' B1{1} '-to-' B1{2}]};
                ind_area.mouse = i;
                ind_area.session= j;
                ind_area.ZT = {['ZT' B2{1} '-to-' B2{2}]};
                tot_mean = ind_mean;
                tot_area = ind_area;
                check=0;
            else
                ind_mean.mouse = i;
                ind_mean.session= j;
                ind_mean.ZT = {['ZT' B1{1} '-to-' B1{2}]};
                ind_area.mouse = i;
                ind_area.session= j;
                ind_area.ZT = {['ZT' B2{1} '-to-' B2{2}]};
                tot_mean(end+1,:) = ind_mean;
                tot_area(end+1,:) = ind_area;
            end
        end
    end
end

%group based on unique ZTs and plot
clf
[group, id] = findgroups(tot_mean.ZT);
gids = unique(group);
for i=1:numel(gids)
    if i==1
        ntot_mean = tot_mean(group==gids(i),:);
        ntot_area = tot_area(group==gids(i),:);
    else
        ntot_mean =[ntot_mean; tot_mean(group==gids(i),:)];
        ntot_area = [ntot_area;tot_area(group==gids(i),:)];
    end
    subplot(2,1,i)
    set(gca,'TickLabelInterpreter', 'none');
    notBoxPlot(table2array(tot_mean(group==gids(i),1:4)),'jitter',0.4);
    xticklabels(tot_mean.Properties.VariableNames(1:4));
    lab_arr = tot_mean.ZT(group==gids(i));
    ylabel(['mean \Deltaf/f ' lab_arr{1}])
end
exportgraphics(gcf,fullfile(dirpath,['mean_seperate_ZTs_all_animals.pdf']),'Resolution',300,'ContentType','vector');
writetable(ntot_mean,fullfile(dirpath,['mean_values_all_animals.csv']));

clf

notBoxPlot(table2array(ntot_mean(:,1:4)),'interval','tinterval')


first_zt = table2array(tot_mean(group==gids(1),1:4));
sec_zt = table2array(tot_mean(group==gids(2),1:4));
first_labels = tot_mean.ZT(group==gids(1));
sec_labels = tot_mean.ZT(group==gids(2));
lab_arr1 = tot_mean.ZT(group==gids(1));
lab_arr2 = tot_mean.ZT(group==gids(2));
clf
%plot pairwise, only works if there are two ZTs'
n=1;
for i=1:4 %there are 4 groups
    plot(ones(length(first_zt(:,i)),1)*n,first_zt(:,i),'.k','MarkerSize',20)
    zt1_coord = [ones(length(first_zt(:,i)),1)*n  first_zt(:,i)];
    zt2_coord = [ones(length(first_zt(:,i)),1)*(n+1) sec_zt(:,i)];
    hold all
    plot(ones(length(first_zt(:,i)),1)*n+1,sec_zt(:,i),'.r','MarkerSize',20)
    plot([zt1_coord(:,1) zt2_coord(:,1)]',[zt1_coord(:,2) zt2_coord(:,2)]','Color',[0.7 0.7 0.7])
    n=n+2;
end
xlim([0 9])
set(gca,'TickLabelInterpreter', 'none');
xticks([1.5:2:8]);
xticklabels(tot_mean.Properties.VariableNames(1:4));
ylabel('mean \Deltaf/f');
legend(lab_arr1{1},lab_arr{2});
exportgraphics(gcf,fullfile(dirpath,['mean_pairwise_ZTs_all_animals.pdf']),'Resolution',300,'ContentType','vector');


%% combine binned means

check = 1;
for i=1:numel(mouse)
    for j=1:numel(mouse(i).session)
        mean_binned_paths = dir(fullfile(mouse(i).session{j},'meanbinnedDFF*'));
        %load all the data
        for ii=1:numel(mean_binned_paths)
            %get ZTs
            B1 = regexp(mean_binned_paths(ii).name,'\d*','Match');
            ind_binned_mean = readtable(fullfile(mean_binned_paths(ii).folder,mean_binned_paths(ii).name));
            if check
                ind_binned_mean.mouse = repmat(i,3,1);
                ind_binned_mean.session= repmat(j,3,1);
                ind_binned_mean.ZT = repmat({['ZT' B1{1} '-to-' B1{2}]},3,1);
                tot_bin_mean = ind_binned_mean;
                check=0;
            else
                ind_binned_mean.mouse = repmat(i,3,1);
                ind_binned_mean.session= repmat(j,3,1);
                ind_binned_mean.ZT = repmat({['ZT' B1{1} '-to-' B1{2}]},3,1);
                tot_bin_mean = [tot_bin_mean ; ind_binned_mean];
            end
        end
    end
end

%group based on unique ZTs and plot
clf
[group, id] = findgroups(tot_bin_mean.ZT);
gids = unique(group);
for i=1:numel(gids)
    if i==1
        ntot_bin_mean = tot_bin_mean(group==gids(i),:);
    else
        ntot_bin_mean =[ntot_bin_mean; tot_bin_mean(group==gids(i),:)];
    end
    subplot(2,1,i)
    set(gca,'TickLabelInterpreter', 'none');
    notBoxPlot(table2array(tot_bin_mean(group==gids(i),1:4)),'jitter',0.4);
    xticklabels(tot_bin_mean.Properties.VariableNames(1:4));
    lab_arr = tot_bin_mean.ZT(group==gids(i));
    ylabel(['mean \Deltaf/f ' lab_arr{1}])
    ylim([-2 5])
end
exportgraphics(gcf,fullfile(dirpath,['mean_binned_seperate_ZTs_all_animals.pdf']),'Resolution',300,'ContentType','vector');
writetable(ntot_bin_mean,fullfile(dirpath,['mean_values_binned_all_animals.csv']));


%% plot transients, can change and re-save as well.

%load previously saved data and reprocess it
bin_reanalyze = 0;
M=define_M();
check = 1;
for i=1:numel(mouse)
    for j=1:numel(mouse(i).session)
        drsfG_Path = dir(fullfile(mouse(i).session{j},'GCaMP_dsrfG*')); %paths to prev processed GCaMP data
        transient_Paths = dir(fullfile(mouse(i).session{j},'all-transients*'));
        %load all the data
        if bin_reanalyze %transient recalculation if needed
            for ii=1:numel(drsfG_Path)
                conv_zt = str2double(regexp(drsfG_Path(ii).name,'\d*','Match'));
                conv_zt = conv_zt(1:2);
                drsfG = readtable(fullfile(drsfG_Path(ii).folder,drsfG_Path(ii).name));
                [allTransients, upTransients] = calc_transients(M,drsfG.dff,drsfG.Score,drsfG.Time,conv_zt,path)
            end
        else
            for ii=1:numel(transient_Paths)
                conv_zt = str2double(regexp(transient_Paths(ii).name,'\d*','Match'));
                ind_transient = readtable(fullfile(transient_Paths(ii).folder,transient_Paths(ii).name));
                %alltransients have Transient Number, Total Secs, Frequency
                if check
                    ind_transient.mouse = repmat(i,3,1);
                    ind_transient.session = repmat(j,3,1);
                    ind_transient.ZT = repmat({['ZT' num2str(conv_zt(1)) '-to-' num2str(conv_zt(2))]},3,1);
                    tot_transient = ind_transient(3,:); %3rd row here is the frequency
                    check=0;
                else
                    ind_transient.mouse = repmat(i,3,1);
                    ind_transient.session = repmat(j,3,1);
                    ind_transient.ZT = repmat({['ZT' num2str(conv_zt(1)) '-to-' num2str(conv_zt(2))]},3,1);
                    tot_transient(end+1,:) = ind_transient(3,:);
                end
                
            end
        end
        
    end
end
clf
[group, id] = findgroups(tot_transient.ZT);
gids = unique(group);
for i=1:numel(gids)
    if i==1
        ntot_transient = tot_transient(group==gids(i),:);
    else
        ntot_transient =[ntot_transient; tot_transient(group==gids(i),:)];
    end
    subplot(2,1,i)
    ylim([0 0.2])
    set(gca,'TickLabelInterpreter', 'none');
    notBoxPlot(table2array(tot_transient(group==gids(i),1:4)),'jitter',0.4);
    xticklabels(tot_transient.Properties.VariableNames(1:4));
    lab_arr = tot_transient.ZT(group==gids(i));
    ylabel(['transient rate (Hz) ' lab_arr{1}])
end
exportgraphics(gcf,fullfile(dirpath,['transients_seperate_ZTs_all_animals.pdf']),'Resolution',300,'ContentType','vector');
writetable(ntot_transient,fullfile(dirpath,['transients_all_animals.csv']));


first_zt = table2array(tot_transient(group==gids(1),1:4));
sec_zt = table2array(tot_transient(group==gids(2),1:4));
first_labels = tot_transient.ZT(group==gids(1));
sec_labels = tot_transient.ZT(group==gids(2));
lab_arr1 = tot_transient.ZT(group==gids(1));
lab_arr2 = tot_transient.ZT(group==gids(2));
clf
%plot pairwise, only works if there are two ZTs'
n=1;
for i=1:4 %there are 4 groups
    plot(ones(length(first_zt(:,i)),1)*n,first_zt(:,i),'.k','MarkerSize',20)
    zt1_coord = [ones(length(first_zt(:,i)),1)*n  first_zt(:,i)];
    zt2_coord = [ones(length(first_zt(:,i)),1)*(n+1) sec_zt(:,i)];
    hold all
    plot(ones(length(first_zt(:,i)),1)*n+1,sec_zt(:,i),'.r','MarkerSize',20)
    plot([zt1_coord(:,1) zt2_coord(:,1)]',[zt1_coord(:,2) zt2_coord(:,2)]','Color',[0.7 0.7 0.7])
    n=n+2;
end
xlim([0 9])
set(gca,'TickLabelInterpreter', 'none');
xticks([1.5:2:8]);
xticklabels(tot_transient.Properties.VariableNames(1:4));
ylabel('transient rate (Hz)');
legend(lab_arr1{1},lab_arr{2});
ylim([0 0.2])
exportgraphics(gcf,fullfile(dirpath,['transients_pairwise_ZTs_all_animals.pdf']),'Resolution',300,'ContentType','vector');

%% calculate all the transitions

%load previously saved transitions
bin_reanalyze = 0;
check = 1;
window = 50;
ds = 4;
for i=1:numel(mouse)
    for j=1:numel(mouse(i).session)
        rsfG_Path = dir(fullfile(mouse(i).session{j},'GCaMP_dff*')); %paths to prev processed GCaMP data
        transition_Paths = dir(fullfile(mouse(i).session{j},'Transition-statistics*'));
        %load all the data
        if bin_reanalyze %transient recalculation if needed
            for ii=1:numel(rsfG_Path)
                rsfG = readtable(fullfile(rsfG_Path(ii).folder,rsfG_Path(ii).name));
                conv_zt = str2double(regexp(rsfG_Path(ii).name,'\d*','Match'));
                conv_zt = conv_zt(1:2);
                path = mouse(i).session{j};
                [~] = calc_transitions(rsfG,ds,window,conv_zt,path);
            end
        else
            for ii=1:numel(transition_Paths)
                clear transition_table %move this into the function
                conv_zt = str2double(regexp(transition_Paths(ii).name,'\d*','Match'));
                load(fullfile(transition_Paths(ii).folder,transition_Paths(ii).name));
                transition_table.traces  = num2cell(m_traces,2);
                transition_table.scores = num2cell(m_scores,2);
                transition_table.tidx = m_tidx';
                transition_table.names = m_name';
                transition_table.numbers = num2cell(m_trans,2);
                transition_table = struct2table(transition_table);
                %alltransients have Transient Number, Total Secs, Frequency
                if check
                    transition_table.mouse = repmat(i,numel(transition_table.tidx),1);
                    transition_table.session = repmat(j,numel(transition_table.tidx),1);
                    transition_table.ZT = repmat({['ZT' num2str(conv_zt(1)) '-to-' num2str(conv_zt(2))]},numel(transition_table.tidx),1);
                    ntransition_table = transition_table; %3rd row here is the frequency
                    check=0;
                else
                    transition_table.mouse = repmat(i,numel(transition_table.tidx),1);
                    transition_table.session = repmat(j,numel(transition_table.tidx),1);
                    transition_table.ZT = repmat({['ZT' num2str(conv_zt(1)) '-to-' num2str(conv_zt(2))]},numel(transition_table.tidx),1);
                    ntransition_table = [ntransition_table ; transition_table];
                end
            end
        end
    end
end

%ntransition_table has all the transitions from all animals
% we don't care about active/quiet wake to NREM transitions,
etransition_table = ntransition_table; %keep the original just in case
%switch all Active to NREM and NREM to Active to Quiet!!
mod_idx = ismember(etransition_table.names,'Active_Wake to NREM');
etransition_table(mod_idx,'names') =  {'Quiet_Wake to NREM'};
etransition_table.numbers(mod_idx) = num2cell([1 2],2);
mod_idx = ismember(etransition_table.names,'NREM to Active_Wake');
etransition_table(mod_idx,'names') =  {'NREM to Quiet_Wake'};
etransition_table.numbers(mod_idx) = num2cell([2 1],2);
%change dscores as well


%desired groups
[group, id] = findgroups(etransition_table.names);

%we are interested in these four groups
desired_groups = [{'Quiet_Wake to NREM'} {'NREM to Quiet_Wake'} {'NREM to REM'} {'REM to NREM'}];


newFolder= fullfile(dirpath,'pca_results');

if exist(newFolder,'dir')
    npath = newFolder;
else
    mkdir(newFolder);
    npath = newFolder;
end
fd =25;
window = 50;
%plot all the transitions 
for i=1:numel(desired_groups)
    sub_t = etransition_table(ismember(etransition_table.names,desired_groups{i}),:);
    %split into ZTs
    [group, id] = findgroups(sub_t.ZT);
    %go through each ZT
    for j=1:numel(gids)
        sub_t1 = sub_t(group==gids(j),:);
        t_arr = cell2mat(sub_t1.traces);
        
        t_arr = smoothdata(t_arr,2,'gaussian',75);
        clf
        subplot(4,1,1:3)
        [coeff,score,latent,tsquared,explained] = pca(t_arr);
        [jnk srtinx] = sort(score(:,1),'descend');
        imagesc(t_arr(srtinx,:))
        tot_sec = size(t_arr,2)/fd;
        xticks([0:fd*10:size(t_arr,2)])
        labelarr = cellfun(@num2str,(num2cell((tot_sec/2)*-1:10:tot_sec/2)),'un',0);
        xticklabels(labelarr)
        ylabel('Transitions')
        xlabel('Time (sec)')
        title(strcat(unique(sub_t1.names),{' '},unique(sub_t1.ZT)),'interpreter','none')
        caxis([-0.5 6])
        subplot(414)
        shadedErrorBar(1:size(t_arr,2),mean(t_arr),std(t_arr,0,1)/sqrt(size(t_arr,1)))
        %make 10 second windows
        xticks([0:fd*10:size(t_arr,2)])
        xticklabels(labelarr)
        xlabel('Time (sec)')
        pname = fullfile(npath,strcat(unique(sub_t1.names),{' '},unique(sub_t1.ZT), {'.pdf'}));
        exportgraphics(gcf,pname{1},'Resolution',300,'ContentType','vector');
        
    end
end



desired_groups = [{'Quiet_Wake to Active_Wake'} {'Active_Wake to Quiet_Wake'}];
fd =25;
window = 50;
%plot all the transitions 
for i=1:numel(desired_groups)
    sub_t = etransition_table(ismember(etransition_table.names,desired_groups{i}),:);
    %split into ZTs
    [group, id] = findgroups(sub_t.ZT);
    %go through each ZT
    for j=1:numel(gids)
        sub_t1 = sub_t(group==gids(j),:);
        t_arr = cell2mat(sub_t1.traces);
        t_arr = smoothdata(t_arr,2,'gaussian',75);
        clf
        subplot(4,1,1:3)
        [coeff,score,latent,tsquared,explained] = pca(t_arr);
        [jnk srtinx] = sort(score(:,1),'descend');
        imagesc(t_arr(srtinx,:))
        tot_sec = size(t_arr,2)/fd;
        xticks([0:fd*10:size(t_arr,2)])
        labelarr = cellfun(@num2str,(num2cell((tot_sec/2)*-1:10:tot_sec/2)),'un',0);
        xticklabels(labelarr)
        ylabel('Transitions')
        xlabel('Time (sec)')
        title(strcat(unique(sub_t1.names),{' '},unique(sub_t1.ZT)),'interpreter','none')
        caxis([0 5])
        subplot(414)
        shadedErrorBar(1:size(t_arr,2),mean(t_arr),std(t_arr,0,1)/sqrt(size(t_arr,1)))
        %make 10 second windows
        xticks([0:fd*10:size(t_arr,2)])
        xticklabels(labelarr)
        xlabel('Time (sec)')
        pname = fullfile(npath,strcat(unique(sub_t1.names),{' '},unique(sub_t1.ZT), {'.pdf'}));
        exportgraphics(gcf,pname{1},'Resolution',300,'ContentType','vector');
        
    end
end

%% calculate binned transients and re-plot it.
%load previously saved data and reprocess it
bin_reanalyze =1;
M=define_M();
check = 1;
bin_plot =1;
for i=1:numel(mouse)
    for j=1:numel(mouse(i).session)
        drsfG_Path = dir(fullfile(mouse(i).session{j},'GCaMP_dsrfG*')); %paths to prev processed GCaMP data
        sT_Path = dir(fullfile(mouse(i).session{j},'EEG-EMG-ZT*'));
        transient_Paths = dir(fullfile(mouse(i).session{j},'all-transients-binned-*'));
        %load all the data
        if bin_reanalyze %transient recalculation if needed
            for ii=1:numel(drsfG_Path)
                conv_zt = str2double(regexp(drsfG_Path(ii).name,'\d*','Match'));
                conv_zt = conv_zt(1:2);
                drsfG = readtable(fullfile(drsfG_Path(ii).folder,drsfG_Path(ii).name));
                sT = readtable(fullfile(sT_Path(ii).folder,sT_Path(ii).name));
                path = mouse(i).session{j};
                [allPeakIds] = calc_transients_v2(M,drsfG.dff,drsfG.Score,drsfG.Time,conv_zt,path);
                if bin_plot
                    [allPeakIds] = plot_windows(M,drsfG.dff,drsfG.Score,drsfG.Time,conv_zt,path);
                end
            end
        else
            for ii=1:numel(transient_Paths)
                conv_zt = str2double(regexp(transient_Paths(ii).name,'\d*','Match'));
                ind_transient = readtable(fullfile(transient_Paths(ii).folder,transient_Paths(ii).name));
                if check
                    ind_transient.mouse = repmat(i,3,1);
                    ind_transient.session = repmat(j,3,1);
                    ind_transient.ZT = repmat({['ZT' num2str(conv_zt(1)) '-to-' num2str(conv_zt(2))]},3,1);
                    tot_transient = ind_transient;
                    check=0;
                else
                    ind_transient.mouse = repmat(i,3,1);
                    ind_transient.session = repmat(j,3,1);
                    ind_transient.ZT = repmat({['ZT' num2str(conv_zt(1)) '-to-' num2str(conv_zt(2))]},3,1);
                    tot_transient = [tot_transient; ind_transient];
                end
                
            end
        end
        
    end
end
clf
[group, id] = findgroups(tot_transient.ZT);
gids = unique(group);
for i=1:numel(gids)
    if i==1
        ntot_transient = tot_transient(group==gids(i),:);
    else
        ntot_transient =[ntot_transient; tot_transient(group==gids(i),:)];
    end
    subplot(2,1,i)
    ylim([0 0.2])
    set(gca,'TickLabelInterpreter', 'none');
    notBoxPlot(table2array(tot_transient(group==gids(i),1:4)),'jitter',0.4);
    xticklabels(tot_transient.Properties.VariableNames(1:4));
    lab_arr = tot_transient.ZT(group==gids(i));
    ylabel(['transient rate (Hz) ' lab_arr{1}])
end
exportgraphics(gcf,fullfile(dirpath,['transients_binned_seperate_ZTs_all_animals.pdf']),'Resolution',300,'ContentType','vector');
writetable(ntot_transient,fullfile(dirpath,['transients_binned_all_animals.csv']));


first_zt = table2array(tot_transient(group==gids(1),1:4));
sec_zt = table2array(tot_transient(group==gids(2),1:4));
first_labels = tot_transient.ZT(group==gids(1));
sec_labels = tot_transient.ZT(group==gids(2));
lab_arr1 = tot_transient.ZT(group==gids(1));
lab_arr2 = tot_transient.ZT(group==gids(2));
clf
%plot pairwise, only works if there are two ZTs'
n=1;
for i=1:4 %there are 4 groups
    plot(ones(length(first_zt(:,i)),1)*n,first_zt(:,i),'.k','MarkerSize',20)
    zt1_coord = [ones(length(first_zt(:,i)),1)*n  first_zt(:,i)];
    zt2_coord = [ones(length(first_zt(:,i)),1)*(n+1) sec_zt(:,i)];
    hold all
    plot(ones(length(first_zt(:,i)),1)*n+1,sec_zt(:,i),'.r','MarkerSize',20)
    plot([zt1_coord(:,1) zt2_coord(:,1)]',[zt1_coord(:,2) zt2_coord(:,2)]','Color',[0.7 0.7 0.7])
    n=n+2;
end
xlim([0 9])
set(gca,'TickLabelInterpreter', 'none');
xticks([1.5:2:8]);
xticklabels(tot_transient.Properties.VariableNames(1:4));
ylabel('transient rate (Hz)');
legend(lab_arr1{1},lab_arr{2});
ylim([0 0.2])
exportgraphics(gcf,fullfile(dirpath,['transients_pairwise_ZTs_all_animals.pdf']),'Resolution',300,'ContentType','vector');





