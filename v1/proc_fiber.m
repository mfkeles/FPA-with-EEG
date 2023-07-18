function chunks = proc_fiber(G)
%processes fiber photometry data that is collected in intervals
time = G.Time_s_(:);

mx_val = max(diff(time));
mx_val = mx_val-10;

slice_idx = find(diff(G.Time_s_) > mx_val);

disp(['Found ' num2str(numel(slice_idx)+1) ,' chunks'])

%add last idx
slice_idx(end+1) = height(G);

for i=1:numel(slice_idx)
    if i==1
        chunks{i} = G(1:slice_idx(i),:);
    else
        chunks{i} = G(slice_idx(i-1)+1:slice_idx(i),:);
    end
end


fs = round( 1/median(diff(chunks{1}.Time_s_))); %calculates the sampling frequency


%make all of the chunks same length - slice 5 seconds from the end
end_point= fs*14*60 + fs*60;

for i=1:numel(chunks)
    chunks{i} = chunks{i}(1:end_point,:);
end


% Correct for movement artifacts.


for i=1:numel(chunks)
    chunks{i}.cF = filt_traces(chunks{i},fs);
end

    function trace = filt_traces(chunk,fs)
        ref_raw  = chunk.AnalogIn__Ch_1AIn_1_Dem_AOut_1_;
        g_raw = chunk.AnalogIn__Ch_1AIn_1_Dem_AOut_2_;

        bleachingFilter = designfilt('lowpassiir', 'HalfPowerFrequency', 0.1, 'SampleRate', fs, 'DesignMethod', 'butter', 'FilterOrder', 12);
        rLowpass = filtfilt(bleachingFilter, ref_raw);
        sLowpass = filtfilt(bleachingFilter, g_raw);
        rFit = fit(chunk.Time_s_, rLowpass, fittype('exp1'));
        sFit = fit(chunk.Time_s_, sLowpass, fittype('exp1'));
        rBleaching = rFit(chunk.Time_s_);
        sBleaching = sFit(chunk.Time_s_);
        rCorrected = ref_raw - rBleaching;
        sCorrected = g_raw - sBleaching;
        trace = sCorrected - rCorrected;
    end
end