function [sT,rsG]= slice_timepoints(time_points,T,G,S,frequency)

start_ind = find(T.Time==time_points{1},1);
stop_ind  = find(T.Time==time_points{2},1);


%EEG+EMG Data Slicing
%find the TTL rise time closest to the given start and stop time.
%get all TTL rise points, each point is a pulse sent by the fiber
%photometry data
[TTL_ind] = find(contains(T.TTL,'TTL: Rise;'));
slice_from = TTL_ind(find(TTL_ind<start_ind,1,'last'));
slice_to = TTL_ind(find(TTL_ind>stop_ind,1));
sT = T(slice_from:slice_to,:);

%Fiber Photometry Data Slicing 
gTTL_ind = find(diff(round(G.DigitalI_O_Ch_1DI_O_1))>0);
gTTL_ind = [1; gTTL_ind]; %we are adding 1 to the first pulse to correctly index the cutoffs from EEG TTL channel
sG = G(gTTL_ind(TTL_ind==slice_from)+1:gTTL_ind(TTL_ind==slice_to)-1,:);


%Using the timestamp column, match the timestamps in the scoring table with
%the timestamps in the EEG+EMG table
[inter_Stamp iEEG iscore] = intersect(sT.TimeStamp,S.TimeStamp);
combined_TimeScores = zeros(size(sT.TimeStamp,1),3);
combined_TimeScores(:,1) = sT.TimeStamp;

for i=1:length(iscore)-1
    epoch_start = find(sT.TimeStamp == S.TimeStamp(iscore(i)),1);
    epoch_end = find(sT.TimeStamp == S.TimeStamp(iscore(i)+1),1);
    scoring = ones(size(epoch_start:epoch_end-1,2),1);
    combined_TimeScores(epoch_start:epoch_end-1,2) = scoring*S.Qiang_0_Numeric(iscore(i));
    combined_TimeScores(epoch_start:epoch_end-1,3) = i;
end

%The beginning and the end of the scores might contain 0s, this is due to
%the fact that EEG+EMG data is 400 Hz and Scoring is 0.2 Hz. So the time
%stamp on the scoring file is not in the EEG+EMG file since we cut that
%based on the TTL. This is not an issue and these score will be ignored
%later. 
sT.Scores = combined_TimeScores(:,2);
sT.Epochs = combined_TimeScores(:,3);


%Epoch indices
e_ind = find(diff(sT.Epochs)==1);
%resample to traget frequency
fs = 1/median(diff(sG.Time_s_)); %sampling frequency
[p ,q] = rat(frequency/fs);
rsG.reference = resample(sG.AnalogIn__Ch_1AIn_1_Dem_AOut_1__LowPass,sG.Time_s_,frequency,p,q);
rsG.signal = resample(sG.AnalogIn__Ch_1AIn_1_Dem_AOut_2__LowPass,sG.Time_s_,frequency,p,q);
rsG.time = resample(sG.AnalogIn__Ch_1AIn_1_Dem_AOut_2__LowPass,sG.Time_s_,frequency,p,q);
rsG.time = resample(sG.Time_s_,sG.Time_s_,frequency,p,q);

%slice the fiber photometry data based on the found indices
fiber_ind = floor((find(diff(sT.Epochs)==1)*length(rsG.time))/length(sT.Epochs));
tmp_arr = [zeros(size(rsG.time,1),1) zeros(size(rsG.time,1),1)]; 
fiber_ind = [0 ; fiber_ind];
e_ind = [0;e_ind];
for i=1:length(fiber_ind)-1
    tmp_arr(fiber_ind(i)+1:fiber_ind(i+1),1) = sT.Epochs(e_ind(i)+1);
    tmp_arr(fiber_ind(i)+1:fiber_ind(i+1),2) = sT.Scores(e_ind(i)+1);
end

rsG.Epochs = tmp_arr(:,1);
rsG.Scores = tmp_arr(:,2);
%remove sections without Epochs
rsG(rsG.Epochs ==0 ,:) =[];
sT(sT.Epochs == 0,:) = [];

end
