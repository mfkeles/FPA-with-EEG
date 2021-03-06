function [M, downsampled, dscore,dtime,depochs]= calc_transitions(rsfG,ds,window,conv_zt,path)

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
    area(k{i}) = (trapz(rsfG.dff(rsfG.Scores == k{i}))/fs);
    headers{k{i}} = M(k{i});
    
end


%find max values for each state binned per hour
k = keys(M);
for i=1:numel(k)
    for j=0:2
        bin_len = floor(length(rsfG.dff)/3);
        bin_rsfG = rsfG((bin_len*j)+1:bin_len*(j+1),:);
        meanbinnedDFF(j+1,k{i}) = mean(bin_rsfG.dff(bin_rsfG.Scores == k{i}));
        areabinned(j+1,k{i}) = (trapz(bin_rsfG.dff(bin_rsfG.Scores == k{i}))/fs);
        headers{k{i}} = M(k{i});
    end
end

meanDFF = array2table(meanDFF);
meanDFF.Properties.VariableNames = headers;
area = array2table(area);
area.Properties.VariableNames = headers;

meanbinnedDFF = array2table(meanbinnedDFF);
meanbinnedDFF.Properties.VariableNames = headers;
meanbinnedDFF.hours = [1:3]';



%downsample the traces by a factor of ds
downsampled = downsample(rsfG.dff,ds);
depochs = downsample(rsfG.Epochs,ds);
dscore = downsample(rsfG.Scores,ds);
dtime = downsample(rsfG.time,ds);

%find all the possible permutations if we pick 2 out of 4 elements
v = perms(keySet);
v = unique(v(:,1:2),'rows');

%all combinations with
for i=1:size(v,1)
    all_combs{i} = [M(v(i,1)) ' to ' M(v(i,2))];
end


%find if there is unscored epochs and make them quietwake by defualt, this
%doesn't affect the meanDFF and other stuff
if find(dscore==255)
    dscore(dscore==255) = 1;
end

%calculate all the transitions
t_idx = find(diff(dscore)~=0); %find all the transitions
for i=1:numel(t_idx)
    transition(i,:) = [dscore(t_idx(i)) dscore(t_idx(i)+1)];
    t_name{i} = [M(transition(i,1)) ' to ' M(transition(i,2))];
end

%composition of transitions
figure
[uni,~,idx] = unique(t_name);
h =histogram(idx,length(unique(idx)));
h.BinEdges = [0.5:1:length(unique(idx))+1];
xticklabels(uni)
set(gca,'TickLabelInterpreter', 'none');
xtickangle(45)
ylabel('Total # of transitions')
no_transition = setdiff(all_combs,t_name);
exportgraphics(gcf,fullfile(path,['Transition_Histogram' '-ZT-' num2str(conv_zt(1)) '-to-' num2str(conv_zt(2)) '.png']),'Resolution',300,'ContentType','image');
close all


skipping = [];
n=1;
%slice all transitions 
for i=1:numel(t_idx)
    if t_idx(i)-window < 0 || t_idx(i)+window > length(downsampled)
        skipping(end+1) = i;
        continue
    else
        m_traces(n,:) = downsampled(t_idx(i)-window+1:(t_idx(i,:)+window));
        m_scores(n,:) = dscore(t_idx(i)-window+1:(t_idx(i,:)+window));
        m_tidx(n) = t_idx(i);
        m_trans(n,:) = transition(i,:);
        m_name{n} = t_name{i};
        n=n+1;
    end
end

drsfG.dff = downsampled;
drsfG.Epochs = depochs;
drsfG.Score = dscore;
drsfG.Time = dtime;
drsfG = struct2table(drsfG);

frequency = fs/ds;
%save the processed files
writetable(drsfG, fullfile(path,['GCaMP_dsrfG' 'ZT-' num2str(conv_zt(1)) '-to-' num2str(conv_zt(2)) '-downsampled-to-' num2str(frequency) 'Hz-' 'processed.csv']));
writetable(meanDFF, fullfile(path,['meanDFF-' 'ZT-' num2str(conv_zt(1)) '-to-' num2str(conv_zt(2)) '.csv']));
writetable(meanbinnedDFF, fullfile(path,['meanbinnedDFF-' 'ZT-' num2str(conv_zt(1)) '-to-' num2str(conv_zt(2)) '.csv']));
save(fullfile(path,['Transition-statistics-' 'ZT-' num2str(conv_zt(1)) '-to-' num2str(conv_zt(2)) '.mat']),'m_traces','m_tidx','m_trans','m_name','no_transition','m_scores')
writetable(area, fullfile(path,['area-' 'ZT-' num2str(conv_zt(1)) '-to-' num2str(conv_zt(2)) '.csv']));
end