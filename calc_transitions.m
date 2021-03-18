function [meanDFF,M, ]= calc_transitions(rsfG,ds,window)

%first define a container that has the keys for 1,2,3 and 4 scoring:
%NonRem                 	    2
%REM                    	    3
%Active Wake                    4
%Quiet Wake                   	1
valueSet = {'REM','NREM','Active_Wake','Quiet_Wake'};
keySet = [3,2,4,1];
M =  containers.Map(keySet,valueSet);
fs = round(1/median(diff(rsfG.time))); %calculates the sampling frequency


%Hz after downsampling since window input is in seconds
window = fs/ds*window;


%find max values for each state
k = keys(M);

for i=1:numel(k)
    meanDFF(k{i}) = mean(rsfG.dff(rsfG.Scores == k{i}));
end

%downsample the traces by a factor of ds
downsampled = downsample(rsfG.dff,ds);
depochs = downsample(rsfG.Epochs,ds);
dscore = downsample(rsfG.Scores,ds);


%find all the possible permutations if we pick 2 out of 4 elements
v = perms(keySet);
v = unique(v(:,1:2),'rows');

%all combinations with
for i=1:size(v,1)
    all_combs{i} = [M(v(i,1)) ' to ' M(v(i,2))];
end

%calculate all the transitions
t_idx = find(diff(dscore)~=0); %find all the transitions
for i=1:numel(t_idx)
    transition(i,:) = [dscore(t_idx(i)) dscore(t_idx(i)+1)];
    t_name{i} = [M(transition(i,1)) ' to ' M(transition(i,2))];
end

%composition of transitions
[uni,~,idx] = unique(t_name);
h =histogram(idx,length(unique(idx)));
h.BinEdges = [0.5:1:length(unique(idx))+1];
xticklabels(uni)
set(gca,'TickLabelInterpreter', 'none');
xtickangle(45)
ylabel('Total # of transitions')
no_transition = setdiff(all_combs,t_name);

skipping = [];
%slice all transitions 
for i=1:numel(t_idx)
    if t_idx(i)-window < 0 || t_idx(i)+window > length(downsampled)
        skipping(end+1) = t_idx(i);
        continue
    else
        t_traces(i,:) = downsampled(t_idx(i)-window+1:(t_idx(i,:)+window));
    end
end


%remove transitions that are not included due to window size
if ~isempty(skipping)
rem_list = t_idx==skipping;
t_idx(rem_list) = [];
transition(rem_list) = [];
t_name(rem_list) = [];
end




if ZT_bin
    save([gdirpath(1:end-4) '_proc_transition_just_wake_ZT12_15.mat'],'comb_dat','comb_scores','comb_depochs')
    disp('saved ZT12-15')
else
    save([gdirpath(1:end-4) '_proc_transition_just_wake_ZT0_3.mat'],'comb_dat','comb_scores','comb_depochs')
    disp('saved ZT0-3')
end

end