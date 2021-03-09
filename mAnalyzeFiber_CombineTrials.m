%point to the folder where the mouse data are located
clear
dirpath = uigetdir;
dirlist = dir(dirpath);
dirlist(1:2)=[];
n=1;
qtoawake =1;
for i=1:length(dirlist)
    if dirlist(i).isdir
        cutdir(n) = dirlist(i);
        n=n+1;
    end
end
mm=1;
for j=1:length(cutdir)
    stop=0;
    %enter folders for each mouse
    mouse_folder = dir(fullfile(cutdir(1).folder,cutdir(j).name));
    mouse_folder(1:2) = [];
    clear zt0_3 zt12_15
    n=1;
    for i=1:length(mouse_folder)
        if mouse_folder(i).isdir
            %find .mat file
            proc_file_zt3 = dir(fullfile(mouse_folder(i).folder,mouse_folder(i).name,'*3.mat'));
            proc_file_zt5 = dir(fullfile(mouse_folder(i).folder,mouse_folder(i).name,'*5.mat'));
            if qtoawake
                for zz = 1:length(proc_file_zt3)
                    if ~isempty(strfind(proc_file_zt3(zz).name,'quiet_vs_active')) && qtoawake
                        proc_file_zt3 = proc_file_zt3(zz);
                        continue
                    end
                end
            end
             if qtoawake
                for zz = 1:length(proc_file_zt5)
                    if ~isempty(strfind(proc_file_zt5(zz).name,'quiet_vs_active')) && qtoawake
                        proc_file_zt5 = proc_file_zt5(zz);
                        continue
                    end
                end
            end
            %load it
            if isempty(proc_file_zt3) || isempty(proc_file_zt5)
                stop = 1;
                break
            end
            load(fullfile(proc_file_zt3.folder,proc_file_zt3.name))
            if n==1
                fn = fieldnames(comb_dat);
                for k=1:numel(fn)
                    zt0_3.(fn{k}) = comb_dat.(fn{k});
                    zt0_3_scores.(fn{k}) = comb_scores.(fn{k});
                    zt0_3_depochs.(fn{k}) = comb_depochs.(fn{k});
                    zt0_3_tryid.(fn{k}) = ones(size(comb_depochs.(fn{k}),1),1);
                end
                n=n+1;
            else
                for k=1:numel(fn)
                    zt0_3.(fn{k}) = [zt0_3.(fn{k}) ;comb_dat.(fn{k})];
                    zt0_3_scores.(fn{k}) = [zt0_3_scores.(fn{k}) ;comb_scores.(fn{k})];
                    zt0_3_depochs.(fn{k}) = [zt0_3_depochs.(fn{k}); comb_depochs.(fn{k})];
                    zt0_3_tryid.(fn{k}) = [zt0_3_tryid.(fn{k}) ;ones(size(comb_depochs.(fn{k}),1),1)*n];
                end
                n=n+1;
            end
        end
    end
    if stop
        continue
    end
    n=1;
    for i=1:length(mouse_folder)
        if mouse_folder(i).isdir
            %find .mat file
            %proc_file_zt5 = dir(fullfile(mouse_folder(i).folder,mouse_folder(i).name,'*5.mat'));
            %load it
            load(fullfile(proc_file_zt5.folder,proc_file_zt5.name))
            if n==1
                fn = fieldnames(comb_dat);
                for k=1:numel(fn)
                    zt12_5.(fn{k}) = comb_dat.(fn{k});
                    zt12_5_scores.(fn{k}) = comb_scores.(fn{k});
                    zt12_5_depochs.(fn{k}) = comb_depochs.(fn{k});
                    zt12_5_tryid.(fn{k}) = ones(size(comb_depochs.(fn{k}),1),1);
                end
                n=n+1;
            else
                for k=1:numel(fn)
                    zt12_5.(fn{k}) = [zt12_5.(fn{k}) ;comb_dat.(fn{k})];
                    zt12_5_scores.(fn{k}) = [zt12_5_scores.(fn{k}) ;comb_scores.(fn{k})];
                    zt12_5_depochs.(fn{k}) = [zt12_5_depochs.(fn{k}); comb_depochs.(fn{k})];
                    zt12_5_tryid.(fn{k}) = [zt12_5_tryid.(fn{k}) ;ones(size(comb_depochs.(fn{k}),1),1)*n];
                end
                n=n+1;
            end
        end
    end
    %create 4 variables to save all the data
    if mm==1
        fn = fieldnames(zt0_3);
        for k=1:numel(fn)
            tot03.(fn{k}) = zt0_3.(fn{k});
            tot1215.(fn{k}) = zt12_5.(fn{k});
            tot03_scores.(fn{k}) = zt0_3_scores.(fn{k});
            tot1215_scores.(fn{k}) = zt12_5_scores.(fn{k});
            tot03_depochs.(fn{k}) = zt0_3_depochs.(fn{k});
            tot1215_depochs.(fn{k}) = zt12_5_depochs.(fn{k});
            tot03_tryid.(fn{k}) = zt0_3_tryid.(fn{k});
            tot1215_tryid.(fn{k}) = zt12_5_tryid.(fn{k});
            tot03_mid.(fn{k}) = ones(size(zt0_3_depochs.(fn{k}),1),1);
            tot1215_mid.(fn{k}) = ones(size(zt0_3_depochs.(fn{k}),1),1);
        end
        mm=mm+1;
    else
        for k=1:numel(fn)
            tot03.(fn{k}) = [tot03.(fn{k})  ; zt0_3.(fn{k})];
            tot1215.(fn{k}) = [tot1215.(fn{k}); zt12_5.(fn{k})];
            tot03_scores.(fn{k}) = [tot03_scores.(fn{k}) ;zt0_3_scores.(fn{k})];
            tot1215_scores.(fn{k}) = [tot1215_scores.(fn{k});zt12_5_scores.(fn{k})];
            tot03_depochs.(fn{k}) = zt0_3_depochs.(fn{k});
            tot1215_depochs.(fn{k}) = zt12_5_depochs.(fn{k});
            tot03_tryid.(fn{k}) = zt0_3_tryid.(fn{k});
            tot1215_tryid.(fn{k}) = zt12_5_tryid.(fn{k});
            tot03_mid.(fn{k}) = [ tot03_mid.(fn{k}) ; ones(size(zt0_3_depochs.(fn{k}),1),1)*mm];
            tot1215_mid.(fn{k}) = [tot1215_mid.(fn{k}) ;ones(size(zt0_3_depochs.(fn{k}),1),1)*mm];
        end
        mm=mm+1;
    end
end

%normalize data 

pause
%process further
clf
for k=1:numel(fn)
    subplot(1,4,k)
    try
    tot_ndat.(fn{k}) = tot03.(fn{k}) - tot03.(fn{k})(:,1);
    plot(mean(tot03.(fn{k})),'LineWidth',1)
    title(fn{k})
    xlim([0 1200])
    xticks([0 600 1200])
    grid on
    ax = gca;
    ax.GridColor = [1, 0, 0];
    ax.GridAlpha = 0.75;
    ax.YGrid = 'off';
    hold all
    catch
        disp('skipping... missing data....')
    end
end

kidx_active2nrem = idx; %first one
kidx_nrem2active = idx; %third cluster
kidx_nrem2rem = idx; %second cluster
kidx_rem2nrem = idx; %second cluster
figure
for k=1:numel(fn)
    subplot(1,4,k)
    imagesc(zt0_3.(fn{k}))
    title(fn{k})
end


clf
clear tot_ndat
for k=1:numel(fn)
    subplot(1,4,k)
    tot_ndat.(fn{k}) = tot1215.(fn{k}) - tot1215.(fn{k})(:,1);
    plot(mean(tot1215.(fn{k})),'LineWidth',1)
    title(fn{k})
    xlim([0 1200])
    xticks([0 600 1200])
    grid on
    ax = gca;
    ax.GridColor = [1, 0, 0];
    ax.GridAlpha = 0.75;
    ax.YGrid = 'off';
    hold all
end

kidx_active2nrem = idx; %second cluster
kidx_nrem2active = idx; %third cluster
kidx_nrem2rem = idx; %1st cluster
kidx_rem2nrem = []; %second cluster
save('zt1215indexes.mat','kidx_active2nrem','kidx_nrem2active','kidx_nrem2rem','kidx_rem2nrem');












