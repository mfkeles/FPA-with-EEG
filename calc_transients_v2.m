function [allPeakIds_binned] = calc_transients_v2(M,downsampled,dscore,dtime,conv_zt,path)
parameters = {2.00, @std, @mean};

%low pass filter for peak detection

fd = round(1/median(diff(dtime)));

%peaksFilter = designfilt('lowpassiir','HalfPowerFrequency',0.2,'SampleRate',fd,'DesignMethod','butter','FilterOrder',12);
 [b a] = butter(2,0.05);
 %peaksSmooth = filtfilt(b,a,downsampled);

%divide it into three min bin

%[~, allPeakIds] = findpeaks(+peaksSmooth,'MinPeakProminence',0.5,'MinPeakWidth',20);

karr = M.keys;
allPeakIds_binned=[];
for ii=0:2
    bin_len = floor(length(downsampled)/3);
    bin_downsampled = downsampled((bin_len*ii)+1:bin_len*(ii+1));
    bin_dscore = dscore((bin_len*ii)+1:bin_len*(ii+1));
    peaksSmooth = filtfilt(b,a,bin_downsampled);
    peakThreshold = threshold(parameters,peaksSmooth);
    if any(peaksSmooth >= peakThreshold)
    [~, allPeakIds] = findpeaks(+peaksSmooth, 'MinPeakHeight', peakThreshold);
    end
   allPeakIds_binned = [allPeakIds_binned ; allPeakIds + (bin_len*ii)];
    for i=1:numel(karr)
        transients(karr{i}) = length(find(bin_dscore(allPeakIds)==karr{i}));
        secs(karr{i}) = round(length(bin_dscore(bin_dscore==karr{i}))/fd);
        freqs(karr{i}) = transients(karr{i})/secs(karr{i});
        all_freqs(ii+1,karr{i}) = freqs(karr{i});
        headers{karr{i}} = M(karr{i});
    end
    
end



allTransients = array2table([all_freqs]);
allTransients.Properties.VariableNames = headers;
allTransients.hours = [1:3]';
writetable(allTransients, fullfile(path,['all-transients-binned-' 'ZT-' num2str(conv_zt(1)) '-to-' num2str(conv_zt(2)) '.csv']));

    function value = threshold(parameters, data)
        % {value1, @mad, @median}
        % {value1, @mad}
        % {value1, @mad, value2}
        % {value1}
        % value
        if iscell(parameters)
            n = numel(parameters);
            if n >= 1
                k = parameters{1};
            else
                k = 2.91;
            end
            if n >= 2
                f2 = parameters{2};
            else
                f2 = @mad;
            end
            if n >= 3
                f3 = parameters{3};
            else
                f3 = @median;
            end
        else
            k = parameters;
            f2 = @mad;
            f3 = @median;
        end
        if isa(f3, 'function_handle')
            value = k * f2(data) + f3(data);
        else
            value = k * f2(data) + f3;
        end
    end
end