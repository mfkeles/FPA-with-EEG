classdef DataProcessor < DataLoader
    properties
        fiber_fs
        proc_data
        score_avgs
        score_summary
        proc_trial_data
        MinPeakDistance = 0.2
        MinPeakProminence = 0.1
        peak_data
    end
    methods
        % Constructor
        function obj = DataProcessor(data_table)
            obj@DataLoader(data_table); % Call the constructor of the parent class DataLoader
        end

        % Function to create a combined EEG and score dataset
        function obj = create_comb_eeg_score(obj)
            % Check if eeg_data and score_data are loaded
            if isempty(obj.eeg_data) || isempty(obj.score_data)
                error('Both EEG and score data must be loaded before calling this function.');
            end

            % Get the number of rows in score_data
            nrows = height(obj.score_data);

            % Initialize the cell array to store expanded score data
            expanded_score = cell(nrows, 1);

            % Iterate through each row of score_data
            for i = 1:nrows
                % Repeat the value in the 5th column 2000 times
                expanded_values = repmat(obj.score_data{i, 5}, 2000, 1);
                % Store the expanded row in the cell array
                expanded_score{i} = expanded_values;
            end

            % Add the cell array with expanded score data as a new column in eeg_data
            obj.eeg_data.ExpandedScore = expanded_score;
        end

        % Function to decimate fiber data before processing
        function obj = decimate_fiber_data(obj, decimation_factor)
            % Check if fiber_data is loaded
            if isempty(obj.fiber_data)
                error('Fiber data must be loaded before calling this function.');
            end

            G = obj.fiber_data;
            GDec = table;

            % Decimate each column in the fiber data using the specified decimation factor
            for i = 1:numel(G.Properties.VariableNames)
                GDec.(G.Properties.VariableNames{i}) = G.(G.Properties.VariableNames{i})(1:decimation_factor:end);
            end

            % Replace the original fiber_data with the decimated fiber data
            obj.fiber_data = GDec;
        end
        function chunks = proc_fiber(obj)
            % Check if fiber_data is loaded
            if isempty(obj.fiber_data)
                error('Fiber data must be loaded before calling this function.');
            end

            G = obj.fiber_data;
            time = G.Time_s_(:);

            mx_val = max(diff(time));
            mx_val = mx_val - 10;

            slice_idx = find(diff(G.Time_s_) > mx_val);

            disp(['Found ' num2str(numel(slice_idx)+1) ,' chunks'])

            % Add last index
            slice_idx(end+1) = height(G);

            for i = 1:numel(slice_idx)
                if i == 1
                    chunks{i} = G(1:slice_idx(i),:);
                else
                    chunks{i} = G(slice_idx(i-1)+1:slice_idx(i),:);
                end
            end

            fs = round(1 / median(diff(chunks{1}.Time_s_))); % Calculate the sampling frequency

            % Make all of the chunks the same length - slice 5 seconds from the end
            end_point = (fs * 14 * 60) + (fs * 60);

            for i = 1:numel(chunks)
                chunks{i} = chunks{i}(1:end_point,:);
            end

            % Correct for movement artifacts
            for i = 1:numel(chunks)
                chunks{i}.cF = filt_traces(chunks{i}, fs);
            end

            % Normalize traces and add to chunks
            for i = 1:numel(chunks)
                for normalizationOption = 1:4
                    normalized_trace = normalize_trace(chunks{i}.cF, chunks{i}.Time_s_, normalizationOption);
                    chunks{i}.(['NormalizedTrace_' num2str(normalizationOption)]) = normalized_trace;
                end
            end

            % Nested function to filter traces
            function trace = filt_traces(chunk, fs)
                ref_raw = chunk.AnalogIn__Ch_1AIn_1_Dem_AOut_1_;
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

            function normalized_trace = normalize_trace(f, time, normalizationOption)
                if normalizationOption == 1
                    % "z-score" ==> (f - mean(f0)) / std(f0)
                    f0 = @mean;
                    f1 = @std;
                elseif normalizationOption == 2
                    % "altered z-score" ==> (f - median(f0)) / median(f0)
                    f0 = @median;
                    f1 = @mad;
                elseif normalizationOption == 3
                    f0 = 0;
                    f1 = 1;
                elseif normalizationOption == 4
                    % Literature:
                    %   df/f ==> (f - f0) / f0
                    %   f0: one of mean, median
                    %   f1: one of std, mad
                    f0 = @median;
                    f1 = @std;
                end

                f0_val = normalize(f0, f, time);
                f1_val = normalize(f1, f, time);
                normalized_trace = (f - f0_val) ./ f1_val;
            end

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


        % Get EEG chunks corresponding to the fiber recordings
        function chunks_table = get_eeg_chunks(obj, chunks, fs_eeg)
            if nargin < 3
                fs_eeg = 400;
            end

            eeg = obj.eeg_data;
            chunks_table = cell(numel(chunks), 1);

            %for each chunk find the corresponding EEG file
            for i = 1:numel(chunks)
                if i == 1
                    stop_time = seconds(chunks{i}.Time_s_(end));
                    [~, ind1] = min(abs(eeg.("Record Time") - stop_time));
                    ind1 = ind1 - 1;

                    %slice and flatten
                    tmptable = table;

                    for m = 1:numel(eeg.Properties.VariableNames)
                        tmpcol = vertcat(eeg.(eeg.Properties.VariableNames{m}){1:ind1});
                        tmptable.(eeg.Properties.VariableNames{m}) = tmpcol;
                    end
                    tmpTime = eeg.("Record Time")(1):seconds(1/fs_eeg):eeg.("Record Time")(ind1 + 1);
                    tmpTime = tmpTime(1:end - 1);
                    tmptable.Time = tmpTime';
                    chunks_table{i} = tmptable;
                else
                    %find start and stop time
                    start_time = seconds(chunks{i}.Time_s_(1));
                    [~, ind1] = min(abs(eeg.("Record Time") - start_time));
                    closest_start = eeg.("Record Time")(ind1);

                    stop_time = seconds(chunks{i}.Time_s_(end));
                    [~, ind2] = min(abs(eeg.("Record Time") - stop_time));
                    ind2 = ind2 - 1;

                    %slice and flatten
                    tmptable = table;

                    for m = 1:numel(eeg.Properties.VariableNames)
                        tmpcol = vertcat(eeg.(eeg.Properties.VariableNames{m}){ind1:ind2});
                        tmptable.(eeg.Properties.VariableNames{m}) = tmpcol;
                    end
                    tmpTime = eeg.("Record Time")(ind1):seconds(1/fs_eeg):eeg.("Record Time")(ind2 + 1);
                    tmpTime = tmpTime(1:end - 1);
                    tmptable.Time = tmpTime';
                    chunks_table{i} = tmptable;
                end
            end
        end

        % Find scores for the fiber data based on the time from start in the score table
        function score_arr = find_scores(obj, fiber_time)
            scores = obj.score_data;
            positions = interp1(scores.TimeFromStart, 1:numel(scores.TimeFromStart), fiber_time);
            score_arr = scores{floor(positions), 5};

            if find(score_arr == 129)
                score_arr(score_arr == 129) = 1;
            end
        end

        % Add Scores column to the fiber_chunks
        function fiber_chunks = add_scores_to_fiber_chunks(obj, fiber_chunks)
            for i = 1:numel(fiber_chunks)
                fiber_chunks{i}.Scores = obj.find_scores(fiber_chunks{i}.Time_s_);
            end
        end


        function obj = process_cell_array(obj,cell_array)
            obj.proc_data = table();
            obj.score_avgs = table();
            obj.proc_trial_data = cell(size(cell_array));  % Initialize new property with the same size as input

            fs = round(1 / median(diff(cell_array{1}.Time_s_)));
            obj.fiber_fs = fs;

            % Loop through each table in the cell array
            for i = 1:length(cell_array)
                current_table = cell_array{i};

                % Discard the first 30 seconds of data
                start_idx = 30 * fs + 1;
                current_table = current_table(start_idx:end, :);

                % Remove negative values from the "NormalizedTrace_2" column
                min_neg_value = min(current_table.NormalizedTrace_2(current_table.NormalizedTrace_2 < 0));
                current_table.NormalizedTrace_2 = current_table.NormalizedTrace_2 - min_neg_value;

                % Add a "Trial" column to the table
                trial_col = repmat(i, height(current_table), 1);
                current_table.Trial = trial_col;

                % Concatenate the processed table with the main processed_data table
                obj.proc_data = [obj.proc_data; current_table(:, {'Time_s_', 'NormalizedTrace_2', 'Scores', 'Trial'})];

                % Store the processed table in the new cell array
                obj.proc_trial_data{i} = current_table;  % Store processed current_table in new property

                % Calculate average value of NormalizedTrace_2 for each score
                score_avg_row = table();

                for score_idx = 1:4
                    if any(current_table.Scores == score_idx)
                        score_rows = current_table.Scores == score_idx;
                        score_avg = mean(current_table.NormalizedTrace_2(score_rows));
                    else
                        score_avg = NaN;
                    end
                    score_avg_row.(sprintf('score_%d', score_idx)) = score_avg;
                end

                % Add a "Trial" column to the score_avg_row
                score_avg_row.Trial = i;

                % Add the score_avg_row to the main score_avgs table
                obj.score_avgs = [obj.score_avgs; score_avg_row];
            end
        end

        function obj = calc_score_summary(obj)
            % Calculate the unique score values and initialize a structure to hold the results
            unique_scores = unique(obj.proc_data.Scores);
            obj.score_summary = struct();

            % Use the provided sampling rate
            fs = obj.fiber_fs;

            % Loop through each unique score value and calculate the total time spent
            for i = 1:length(unique_scores)
                current_score = unique_scores(i);
                score_rows = obj.proc_data.Scores == current_score;
                time_spent = sum(score_rows) / fs / 60; % Convert to minutes
                obj.score_summary.(sprintf('score_%d', current_score)) = time_spent;
            end
        end
        function result = calculatePeaks(obj)
            % Initialize the result table
            result = table();

            % Loop through each table in the cell array
            for i = 1:length(obj.proc_trial_data)
                current_table = obj.proc_trial_data{i};

                % Find peaks in the NormalizedTrace_2 data
                [pks,locs,w,p] = findpeaks(current_table.NormalizedTrace_2, 'MinPeakDistance', obj.fiber_fs * obj.MinPeakDistance, 'MinPeakProminence', obj.MinPeakProminence);

                % Initialize arrays to store results for each score
                peak_counts = zeros(1,4);
                peak_avgs = zeros(1,4);
                width_avgs = zeros(1,4);
                prom_avgs = zeros(1,4);

                % Loop through each score and calculate statistics
                for score_idx = 1:4
                    if any(current_table.Scores == score_idx)
                        score_rows = current_table.Scores == score_idx;
                        score_duration = height(current_table(score_rows,:)) / obj.fiber_fs;
                        score_rows_indices = find(score_rows);
                        score_peaks = ismember(locs, score_rows_indices);

                        % Calculate frequency of peaks (number of peaks/duration)
                        peak_counts(score_idx) = sum(score_peaks) / score_duration;

                        % Calculate average peak value, width and prominence
                        peak_avgs(score_idx) = mean(pks(score_peaks));
                        width_avgs(score_idx) = mean(w(score_peaks));
                        prom_avgs(score_idx) = mean(p(score_peaks));
                    else
                        peak_counts(score_idx) = NaN;
                        peak_avgs(score_idx) = NaN;
                        width_avgs(score_idx) = NaN;
                        prom_avgs(score_idx) = NaN;
                    end
                end

                % Create a row for the current trial
                trial_row = array2table([peak_counts, peak_avgs, width_avgs, prom_avgs], 'VariableNames', ...
                    [strcat('peak_count_score_', string(1:4)), strcat('peak_avg_score_', string(1:4)), ...
                    strcat('width_avg_score_', string(1:4)), strcat('prom_avg_score_', string(1:4))]);
                trial_row.Trial = i;

                % Add the row to the result table
                result = [result; trial_row];
            end
        end
        function transitionTable = extractTransitions(obj, persistence)
            % Convert persistence from seconds to number of samples
            persistence = round(persistence * obj.fiber_fs);

            % Initialize empty table for transitions
            transitionTable = table();

            % Iterate over each table in the cell array
            for i = 1:length(obj.proc_trial_data)
                data = obj.proc_trial_data{i};
                scores = data.Scores;
                normalizedTrace = data.NormalizedTrace_2;

                % Initialize variables to track the current state and its start time
                currentState = scores(1);
                stateStartIndex = 1;
                for j = 2:height(data)
                    % Check if the state has changed
                    if scores(j) ~= currentState
                        % Check if the previous state persisted long enough
                        if j - stateStartIndex >= persistence
                            % Check if the next state persists long enough
                            nextStateEndIndex = j + find(scores(j+1:end) ~= scores(j), 1) - 1;
                            if isempty(nextStateEndIndex)
                                nextStateEndIndex = height(data);
                            end
                            if nextStateEndIndex - j >= persistence
                                % If both states persisted long enough, add the corresponding portion of the normalized trace to the transition table
                                transition = normalizedTrace(j-persistence:j+persistence-1);
                                transitionScores = scores(j-persistence:j+persistence-1);
                                transitionType = sprintf('%d to %d', currentState, scores(j));
                                transitionRow = table({transition}, {transitionScores}, {transitionType}, 'VariableNames', {'Transition', 'TransitionScores', 'TransitionType'});
                                transitionTable = [transitionTable; transitionRow];
                            end
                        end
                        % Update the current state and its start time
                        currentState = scores(j);
                        stateStartIndex = j;
                    end
                end
            end
        end    
    end
end

