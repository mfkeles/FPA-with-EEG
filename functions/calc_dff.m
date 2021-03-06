function [rsfG] = calc_dff(rsG,normalizationOption)


ref_raw = rsG.reference;
g_raw = rsG.signal;
fs = round( 1/median(diff(rsG.time))); %calculates the sampling frequency

bleachingFilter = designfilt('lowpassiir', 'HalfPowerFrequency', 0.1, 'SampleRate', fs, 'DesignMethod', 'butter', 'FilterOrder', 12);
rLowpass = filtfilt(bleachingFilter, ref_raw);
sLowpass = filtfilt(bleachingFilter, g_raw);
rFit = fit(rsG.time, rLowpass, fittype('exp1'));
sFit = fit(rsG.time, sLowpass, fittype('exp1'));
rBleaching = rFit(rsG.time);
sBleaching = sFit(rsG.time);
rCorrected = ref_raw - rBleaching;
sCorrected = g_raw - sBleaching;

% Correct for movement artifacts.
f = sCorrected - rCorrected;

lowpassFilter = designfilt('lowpassiir', 'HalfPowerFrequency', 2, 'SampleRate', fs, 'DesignMethod', 'butter', 'FilterOrder', 12);
fSmooth = filtfilt(lowpassFilter, f);

if normalizationOption==1
    % "z-score" ==> (f - mean(f0)) / std(f0)
    f0 = @mean;
    f1 = @std;
elseif  normalizationOption==2
    % "altered z-score" ==> (f - median(f0)) / median(f0)
    f0 = @median;
    f1 = @mad;
elseif normalizationOption==3
    % Doric_photom_analysis.m if plotting variable "normDat".
    f0 = 0;
    f1 = 1;
elseif normalizationOption==4
    % Literature:
    %   df/f ==> (f - f0) / f0
    %   f0: one of mean, median
    %   f1: one of std, mad
    f0 = @median;
    f1 = @std;
end

f0 = normalize(f0,fSmooth,rsG.time);
f1 = normalize(f1,fSmooth,rsG.time);
rsG.dff = (fSmooth -f0)./f1;
rsfG = rsG;
    function output = normalize(parameters, f, time)
        if iscell(parameters)
            fcn = parameters{1};
            if numel(parameters) == 1
                parameters{2} = [-Inf, Inf];
            end
            if isscalar(parameters{2})
                % Produce a vector from moving window.
                if numel(parameters) <= 2.
                    options = {'EndPoints', 'shrink'};
                else
                    options = parameters(3:end);
                end
                frequency = 1 / median(diff(time));
                nSamples = numel(time);
                window = parameters{2};
                window = min(round(window * frequency), nSamples);
                output = fcn(f, window, options{:});
            else
                % Produce a value from all data (or epochs).
                epochs = parameters{2};
                ids = time2id(time, epochs);
                output = fcn(f(ids));
            end
        elseif isa(parameters, 'function_handle')
            % Produce a value from all data (or epochs).
            fcn = parameters;
            epochs = [-Inf, Inf];
            ids = time2id(time, epochs);
            output = fcn(f(ids));
        else
            output = parameters;
        end
    end
end