function [uniquePeakIds,allPeak] = calc_transients(M,downsampled,dscore,dtime)

parameters = {2.00, @std, @mean};

%low pass filter for peak detection

fd = round(1/median(diff(dtime)));

peaksFilter = designfilt('lowpassiir','HalfPowerFrequency',0.2,'SampleRate',fd,'DesignMethod','butter','FilterOrder',12);
peaksSmooth = filtfilt(peaksFilter,downsampled);

peakThreshold = threshold(parameters,peaksSmooth);

if any(peaksSmooth >= peakThreshold)
    [~, uniquePeakIds] = findpeaks(+peaksSmooth, 'MinPeakHeight', peakThreshold);
end

[~, allPeakIds] = findpeaks(+peaksSmooth,'MinPeakProminence',0.5,'MinPeakWidth',20);

karr = M.keys;
for i=1:numel(karr)
    transients(karr{i}) = length(find(dscore(allPeakIds)==karr{i}));
    secs(karr{i}) = round(length(dscore(dscore==karr{i}))/fd);
    freqs(karr{i}) = transients(karr{i})/secs(karr{i});
    headers{karr{i}} = M(karr{i});
end

allTransients = array2table([transients;secs;freqs]);
allTransients.Properties.VariableNames = headers;
allTransients.Row = {'Transient Number','Total Secs','Frequency'};


karr = M.keys;
for i=1:numel(karr)
    utransients(karr{i}) = length(find(dscore(uniquePeakIds)==karr{i}));
    usecs(karr{i}) = round(length(dscore(dscore==karr{i}))/fd);
    ufreqs(karr{i}) = utransients(karr{i})/usecs(karr{i});
    headers{karr{i}} = M(karr{i});
end
upTransients = array2table([utransients;usecs;ufreqs]);
upTransients.Properties.VariableNames = headers;
upTransients.Row = {'Transient Number','Total Secs','Frequency'};



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