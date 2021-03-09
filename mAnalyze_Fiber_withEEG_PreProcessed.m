clear
[gfile, path] = uigetfile('*.csv','Pick .csv file for GCAMP signal');
cd(path)
%check if sliced files exist
sliced_path = dir(fullfile(path,'*sliced.txt'));

if ~isempty(sliced_path)
    disp('Previously processed files found, will load them...')
    for i=1:length(sliced_path)
        if sliced_path(i).name(1:3) == 'EEG'
            if ~isempty(strfind(sliced_path(i).name,'ZT0'))
                ZT0_3_EEG = sliced_path(i);
            else
                ZT12_15_EEG = sliced_path(i);
            end
        else
            if ~isempty(strfind(sliced_path(i).name,'ZT0'))
                ZT0_3_gcamp = sliced_path(i);
            else
                ZT12_15_gcamp = sliced_path(i);
            end
        end
    end
else
    disp('No previously processed files found...')
    gdirpath = fullfile(path,gfile);
    cd(path)
    [efile, path] = uigetfile('*.txt', 'Pick your EEG file');
    edirpath = fullfile(path,efile);
    [efile, path] = uigetfile('*.edf', 'Pick your EEG.edf file');
    edirpath = fullfile(path,efile);
    [sfile, path] = uigetfile('*.txt', 'Pick your score file');
    sdirpath = fullfile(path,sfile);
end


paths = [ZT0_3_gcamp ZT12_15_gcamp];

for n=1:length(paths)
    %load ZT0_3 EEG and GCaMP data
    clear comb_ind comb_dat comb_scores comb_epochs
    modFiber = readtable(fullfile(paths(n).folder,paths(n).name));
    downsampled = downsample(modFiber.dffLowpass,50);
    depochs = downsample(modFiber.Epochs,50);
    dscore = downsample(modFiber.Scores,50);
    %NonRem                 	2
    %REM                    	3
    %Active Wake              	4
    %Quiet Wake                   	1
    %find where quiet wakes are
    [row,col] = find(dscore==1);
    [drow ,dcol] = find(diff(row)>1);
    %drow are the locations of NREM switches.
    %check what comes before it to seperate
    clear comb_ind
    comb_ind.qwake2rem = [];
    comb_ind.rem2qwake = [];
    comb_ind.qwake2awake = [];
    comb_ind.awake2qwake = [];
    for i=1:length(drow)
        %after nrem
        if dscore(row(drow(i))-5) == 1 %leaving quiet wake
            if dscore(row(drow(i))+5) == 3  %quiet wake to REM
                comb_ind.qwake2rem(end+1) = row(drow(i));
            elseif dscore(row(drow(i))+5) ==4 %
                comb_ind.qwake2awake(end+1) = row(drow(i));
            end
        end
        %before nrem
        if dscore(row(drow(i)+1)-5)==  3 %REM to awake
            if dscore(row(drow(i)+1)+5) == 1 %internal check that it's quiet wake
                comb_ind.rem2qwake(end+1) = row(drow(i)+1);
            end
        elseif dscore(row(drow(i)+1)-5) ==4 %active to quiet
            if dscore(row(drow(i)+1)+5) == 1 %internal check
                comb_ind.awake2qwake(end+1) = row(drow(i)+1);
            end
        end
    end
    %add +/- 24 seconds
    fn = fieldnames(comb_ind);
    clear comb_dat comb_scores comb_depochs
    for k=1:numel(fn)
        comb_dat.(fn{k}) = [];
        comb_scores.(fn{k}) = [];
        comb_depochs.(fn{k}) = [];
        %comb_dat(k).scores = [];
        %comb_dat(k).ident = fn{k};
        for i=1:length(comb_ind.(fn{k}))
            if comb_ind.(fn{k})(i)-600<0 || comb_ind.(fn{k})(i)+600 > length(downsampled)
                continue
            else
                comb_dat.(fn{k})(end+1,:) = downsampled(comb_ind.(fn{k})(i)-600:comb_ind.(fn{k})(i)+600);
                comb_scores.(fn{k})(end+1,:) = dscore(comb_ind.(fn{k})(i)-600:comb_ind.(fn{k})(i)+600);
                comb_depochs.(fn{k})(end+1,:) = depochs(comb_ind.(fn{k})(i)-600:comb_ind.(fn{k})(i)+600);
            end
        end
    end
    if isempty(strfind(paths(n).name,'ZT0_3'))
        save(['proc_transition_quiet_vs_active_ZT12_15.mat'],'comb_dat','comb_scores','comb_depochs')
        disp('saved ZT12-15')
    else
        save(['proc_transition_quiet_vs_active_ZT0_3.mat'],'comb_dat','comb_scores','comb_depochs')
        disp('saved ZT0-3')
    end
    
end

