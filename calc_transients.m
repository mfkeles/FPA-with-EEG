function comb_dat = calc_transients(modFiber)

%first define a container that has the keys for 1,2,3 and 4 scoring:
%NonRem                 	    2
%REM                    	    3
%Active Wake                    4
%Quiet Wake                   	1
keySet = {'REM','NREM','Active_Wake','Quiet_Wake'};
valueSet = [3,2,4,1];

M =  containers.Map(keySet,valueSet);





%downsample it by a factor of 50, brings it down to 1250/50 to 25Hz
downsampled = downsample(modFiber.dffLowpass,50);
depochs = downsample(modFiber.Epochs,50);
dscore = downsample(modFiber.Scores,50);

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