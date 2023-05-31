classdef Plotter < handle
    properties
        cmap = [102 194 165
            252 141 98
            141 160 203
            231 138 195] / 255;
    end
    methods
        function obj = Plotter()
            % Constructor
        end
        function plot_fiber_chunks_with_scores(obj, row, fiber_chunks_with_scores, fiber_dir)
            % Get the parent folder of the fiber file
            parent_folder = fileparts(fiber_dir);

            % Create the figure
            fig = figure;

            % Set figure size to A4 dimensions
            set(fig, 'Units', 'centimeters', 'PaperUnits', 'centimeters');
            set(fig, 'PaperSize', [21, 29.7]); % A4 dimensions in centimeters
            set(fig, 'PaperPosition', [0, 0, 21, 29.7]);

            % Plot the data
            for i = 1:numel(fiber_chunks_with_scores)
                subplot(numel(fiber_chunks_with_scores), 1, i);
                chunk = fiber_chunks_with_scores{i};
                fs_fiber = round(1 / median(diff(chunk.Time_s_)));
                start_idx = fs_fiber * 2 + 1; % Skip the first 2 seconds

                patch([start_idx:length(chunk.cF) nan], [chunk.cF(start_idx:end)' nan], [chunk.Scores(start_idx:end)' nan], [chunk.Scores(start_idx:end)' nan], 'edgecolor', 'interp');

                % Set x-axis ticks
                xticks(0:fs_fiber*60:length(chunk.cF));
                xticklabels(0:1:floor(length(chunk.cF)/(fs_fiber*60)));
                xlabel('Time (minutes)');

                % Set title for each subplot
                title(sprintf('%s %s: Trial %d', row.mouse_dir{1}, row.trial_dir{1}, i));
                caxis([1 4])
            end

            colormap(obj.cmap)
            caxis([1 4])

            % Create a legend
            hold on;
            h = zeros(4, 1);
            for i = 1:4
                h(i) = plot(NaN, NaN, 'color', obj.cmap(i, :), 'LineWidth', 2);
            end
            legend(h, 'Wake', 'NREM', 'REM', 'Awake');
            hold off;

            % Save the figure to the 'figs' subfolder
            if ~exist(fullfile(parent_folder, 'figs'), 'dir')
                mkdir(fullfile(parent_folder, 'figs'));
            end
            saveas(fig, fullfile(parent_folder, 'figs', 'fiber_chunks.png'));
            close(fig);
        end

        function plot_summary(obj, data_processor, parent_folder)
            % Get the score_summary from data_processor
            score_summary = data_processor.score_summary;

            % Plot the summary figure
            fig = figure;
            score_values = [score_summary.score_1, score_summary.score_2, score_summary.score_3, score_summary.score_4];
            unique_scores = 1:4;
            bar_handle = bar(unique_scores, score_values, 'FaceColor', 'flat');

            % Set colors for each bar individually
            for i = 1:length(unique_scores)
                bar_handle.CData(i, :) = obj.cmap(i, :);
            end

            % Set x-axis labels
            xticks(unique_scores);
            xticklabels({'Wake', 'NREM', 'REM', 'Awake'});

            % Set y-axis label
            ylabel('Time Spent (minutes)');

            % Set title
            title('Summary of Time Spent in Each Score Value');

            % Save the figure to the 'figs' subfolder
            if ~exist(fullfile(parent_folder, 'figs'), 'dir')
                mkdir(fullfile(parent_folder, 'figs'));
            end
            saveas(fig, fullfile(parent_folder, 'figs', 'summary.png'));
            close(fig);
        end
        function plot_avg_z_scores_per_mouse(obj, combined_data, parent_folder)
            % Set the folder path for saving the plot
            folder_path = 'X:\0-mWAKE-Cre-Gcamp-EEG-NEW';

            % Find unique mouse IDs
            unique_mice = unique(combined_data.mouse);

            % Initialize a table to store average values per mouse
            avg_per_mouse = array2table(zeros(length(unique_mice), 4), 'VariableNames', {'score_1', 'score_2', 'score_3', 'score_4'});

            % Calculate average values for each mouse
            for i = 1:length(unique_mice)
                mouse_id = unique_mice(i);
                mouse_data = combined_data(combined_data.mouse == mouse_id, :);

                % Calculate the mean for each score column, ignoring NaN values
                avg_per_mouse{i, 'score_1'} = nanmean(mouse_data.score_1);
                avg_per_mouse{i, 'score_2'} = nanmean(mouse_data.score_2);
                avg_per_mouse{i, 'score_3'} = nanmean(mouse_data.score_3);
                avg_per_mouse{i, 'score_4'} = nanmean(mouse_data.score_4);
            end

            % Plot the average values per mouse
            fig = figure;
            bar_handle = bar(avg_per_mouse.Variables, 'grouped');

            % Set colors for each bar individually
            for i = 1:size(obj.cmap, 1)
                bar_handle(i).FaceColor = obj.cmap(i, :);
            end

            % Set x-axis labels
            xticks(1:length(unique_mice));
            xticklabels(arrayfun(@(x) sprintf('Mouse %d', x), unique_mice, 'UniformOutput', false));

            % Set y-axis label
            ylabel('Z-score');

            % Set title
            title('Average Z-scores per Mouse');

            % Add legend
            legend({'Score 1 (Wake)', 'Score 2 (NREM)', 'Score 3 (REM)', 'Score 4 (Awake)'}, 'Location', 'northoutside');

            % Save the figure to the specified folder path
            saveas(fig, fullfile(folder_path, 'average_z_scores_per_mouse.png'));
            close(fig);
        end

        function plotZTScores(obj, combined_data, parent_folder)
            % Get unique mouse IDs
            unique_mice = unique(combined_data.mouse);

            % Calculate the mean values for ZT3 and ZT15 separately
            zt3_rows = combined_data.ZT == "3";
            zt15_rows = combined_data.ZT == "15" | combined_data.ZT=="14";

            % Initialize the mean values array
            mean_values_zt3 = zeros(length(unique_mice), 4);
            mean_values_zt15 = zeros(length(unique_mice), 4);

            % Loop through each unique mouse ID and calculate the mean values
            for i = 1:length(unique_mice)
                mouse_id = unique_mice(i);
                mouse_rows = combined_data.mouse == mouse_id;

                % Calculate the mean values for ZT3 and ZT15 for the current mouse
                mean_values_zt3(i, :) = nanmean(combined_data{mouse_rows & zt3_rows, {'score_1', 'score_2', 'score_3', 'score_4'}});
                mean_values_zt15(i, :) = nanmean(combined_data{mouse_rows & zt15_rows, {'score_1', 'score_2', 'score_3', 'score_4'}});
            end

            % Set up the figure
            fig = figure;

            % Set colors for each score type
            cmap = obj.cmap;

            % Define the markers for ZT3 and ZT14-15
            markers = 'o+';

            % Create dummy plots for legend
            for score_idx = 1:4
                scatter(NaN, NaN, 50, 'o', 'MarkerFaceColor', cmap(score_idx, :), 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
                hold on;
                scatter(NaN, NaN, 50, '+', 'MarkerEdgeColor', cmap(score_idx, :), 'LineWidth', 1.5);
            end

            % Loop over each mouse
            for mouse_idx = 1:length(unique_mice)
                % Plot ZT3 and ZT14-15 values for each score
                for score_idx = 1:4
                    x = (mouse_idx-1)*5 + score_idx;  % Adjust x for each score type to avoid overlap
                    y = [mean_values_zt3(mouse_idx, score_idx), mean_values_zt15(mouse_idx, score_idx)];  % ZT3 and ZT14-15 values
                    % ZT3 and ZT14-15 values

                    % Plot ZT3 value
                    scatter(x, y(1), 50, 'o', 'MarkerFaceColor', cmap(score_idx, :), 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
                    hold on;

                    % Plot ZT14-15 value
                    scatter(x, y(2), 50, '+', 'MarkerEdgeColor', cmap(score_idx, :), 'LineWidth', 1.5);

                    % Plot a line connecting the ZT3 and ZT14-15 values
                    plot(x*ones(size(y)), y, 'Color', cmap(score_idx, :), 'LineWidth', 1.5);
                end
            end

            % Adjust x-axis
            xticks(1:5:length(unique_mice)*5);
            xticklabels(arrayfun(@(x) sprintf('Mouse %d', x), unique_mice, 'UniformOutput', false));

            % Set y-axis label
            ylabel('Z-score');

            % Set title
            title('ZT3 and ZT14-15 Z-scores per Mouse');

            % Add legend
            legend({'Score 1 (Wake) - ZT3', 'Score 1 (Wake) - ZT14-15', 'Score 2 (NREM) - ZT3', 'Score 2 (NREM) - ZT14-15', 'Score 3 (REM) - ZT3', 'Score 3 (REM) - ZT14-15', 'Score 4 (Awake) - ZT3', 'Score 4 (Awake) - ZT14-15'}, 'Location', 'northoutside', 'NumColumns', 2);

            % Save the figure
            saveas(fig, fullfile(parent_folder, 'zt3_zt15_zscores_per_mouse.png'));
            close(fig);

        end

        function plotDataChunks(obj, fiber_chunks_with_scores, eeg_chunks_table, fiber_sampling_frequency, eeg_sampling_frequency, parent_folder, decimation_factor)
            fs_eeg = eeg_sampling_frequency;
            MinPeakDistance = 0.2;
            MinPeakProminence = 0.4;
            MinPeakHeight = 0.2;

            % Check if both cell arrays are of the same length
            assert(length(fiber_chunks_with_scores) == length(eeg_chunks_table), ...
                'fiber_chunks_with_scores and eeg_chunks_table should be of same length');

            % Loop over each cell of the cell arrays
            for i = 1:length(fiber_chunks_with_scores)

                fiber_table = fiber_chunks_with_scores{i};
                eeg_table = eeg_chunks_table{i};

                % Remove the first and last 30 seconds
                fiber_table = fiber_table(30*fiber_sampling_frequency+1:end-30*fiber_sampling_frequency, :);
                eeg_table = eeg_table(30*eeg_sampling_frequency+1:end-30*eeg_sampling_frequency, :);

                % Calculate total duration in minutes
                fiber_duration = height(fiber_table) / (60 * fiber_sampling_frequency);  % duration in minutes for fiber data
                eeg_duration = height(eeg_table) / (60 * eeg_sampling_frequency);  % duration in minutes for eeg data

                % Create a new figure for each cell
                %fig = figure('PaperSize',[8.3, 14.7]); % A4 paper size

                % Loop over each bin
                for j = 1:2  % Adjust the loop to iterate over 2 chunks
                    % Calculate start and end indices for the current bin
                    fiber_start_index = round(1 + (j-1) * fiber_duration/2 * 60 * fiber_sampling_frequency);  % Adjust to divide by 2 instead of 3
                    fiber_end_index = round(j * fiber_duration/2 * 60 * fiber_sampling_frequency);  % Adjust to divide by 2 instead of 3
                    eeg_start_index = round(1 + (j-1) * eeg_duration/2 * 60 * eeg_sampling_frequency);  % Adjust to divide by 2 instead of 3
                    eeg_end_index = round(j * eeg_duration/2 * 60 * eeg_sampling_frequency);  % Adjust to divide by 2 instead of 3

                    % Extract the data for the current bin
                    fiber_bin = fiber_table(fiber_start_index:fiber_end_index, :);
                    eeg_bin = eeg_table(eeg_start_index:eeg_end_index, :);

                    % Find peaks for fiber_bin
                    dec_fiber_bin = decimate(fiber_bin.NormalizedTrace_2, decimation_factor);
                    dec_scores = decimate(fiber_bin.Scores, decimation_factor);

                    dec_EMG = decimate(eeg_bin.EMGEMG, decimation_factor);
                    dec_EMG_scores = decimate(eeg_bin.ExpandedScore, decimation_factor);

                    [pks,locs] = findpeaks(dec_fiber_bin- (min(dec_fiber_bin)), 'MinPeakDistance', fiber_sampling_frequency/decimation_factor * MinPeakDistance, 'MinPeakProminence', MinPeakProminence,'WidthReference','halfheight','MinPeakHeight',MinPeakHeight);
                    % Create subplots
                    subplot(6, 1, j*3-2) % first row for fiber data
                    patch([1:length(dec_fiber_bin), nan], [dec_fiber_bin', nan], ...
                        [dec_scores', nan], 'EdgeColor', 'interp','FaceColor','none','LineWidth',1);

                    colormap(gca, obj.cmap);  % Use the colormap defined in Plotter class for fiber data
                    caxis(gca, [1 4]);  % Set color axis limits to match score values
                    hold on;

                    % Then plot the peaks using 'v' as the marker (an upside down triangle)
                    plot(locs, dec_fiber_bin(locs), 'o','Color','r','MarkerSize',3);
                    hold off;
                    title(['Fiber Data Chunk ' num2str(j)]);
                    xticks(0:fiber_sampling_frequency/decimation_factor*60:length(dec_fiber_bin));
                    xticklabels(0:length(dec_fiber_bin)/fiber_sampling_frequency/decimation_factor);
                    xlim([0,length(dec_fiber_bin)])
                    xlabel('Time (min)');
                    ylabel('Normalized Trace');

                    subplot(6, 1, j*3-1); % second row for EEG data
                    pspectrum(eeg_bin.EEGEEG2A_B,fs_eeg,'spectrogram','FrequencyLimits',[0 25],'Leakage',0.85,'OverlapPercent',60,'TimeResolution',10);
                    % After your call to pspectrum
                    colorbarHandles = findobj(gcf,'type','colorbar');
                    delete(colorbarHandles);
                    colormap(gca, 'default');  % Use the default colormap for EEG data
                    title(['EEG Data Chunk ' num2str(j)]);
                    %xlim([0,length(eeg_bin.EEGEEG2A_B)])
                    xlabel('Time (min)');
                    ylabel('EEG Data');

                    subplot(6, 1, j*3) % third row for EMG data
                    patch([1:length(eeg_bin.EMGEMG), nan], [eeg_bin.EMGEMG', nan], ...
                        [eeg_bin.ExpandedScore', nan], 'EdgeColor', 'interp','FaceColor','none','LineWidth',1);

                    colormap(gca, obj.cmap);  % Use the colormap defined in Plotter class for EMG data
                    caxis(gca, [1 4]);  % Set color axis limits to match score values
                    xticks(0:fs_eeg*60:length(eeg_bin.EMGEMG));
                    xticklabels(0:length(eeg_bin.EMGEMG)/fs_eeg);
                    xlim([0,length(eeg_bin.EMGEMG)])
                    ylim([-500,500])
                    title(['EMG Data Chunk ' num2str(j)]);
                    xlabel('Time (min)');
                    ylabel('EMG Data');
                end

                % Save the figure to the 'figs' subfolder
                if ~exist(fullfile(parent_folder, 'figs'), 'dir')
                    mkdir(fullfile(parent_folder, 'figs'));
                end
                %plot2svg(fullfile(parent_folder, 'figs', strcat('fiber_chunks_combined_eeg_trial_pn_',num2str(i),'_patch.svg')))
                exportgraphics(gcf, fullfile(parent_folder, 'figs', strcat('fiber_chunks_combined_eeg_trial_pn_',num2str(i),'_patch_cbar.pdf')),'ContentType','vector','Resolution',150);
                %close(fig);
                clf
            end
        end



    end
    methods (Static)
        function plotTransitionHeatmaps(combined_transitions, combineWakeAwake,norm,fs_fiber,parent_folder)
            % Create a map for scores to states
            scoreToState = {'Wake', 'NREM', 'REM', 'Awake'};

            % If combineWakeAwake is true, combine '1' and '4' into 'cwake'
            if combineWakeAwake
                combined_transitions.TransitionType = regexprep(combined_transitions.TransitionType, '1|4', 'cwake');
                scoreToState = {'cwake', 'NREM', 'REM'};
            end

            % Extract unique transitions
            unique_transitions = unique(combined_transitions.TransitionType);

            % Loop over each unique transition
            for i = 1:length(unique_transitions)
                transition_type = unique_transitions(i);

                % Filter rows for current transition type
                transition_rows = strcmp(combined_transitions.TransitionType, transition_type);

                % Get transition data for current transition type
                transition_data = combined_transitions.Transition(transition_rows);

                % Convert cell array to matrix for plotting
                transition_matrix = horzcat(transition_data{:});
                transition_matrix = transition_matrix';


                if norm
                    % Normalize each row so that the max of each row is equal to 1
                    transition_matrix = transition_matrix ./ max(transition_matrix, [], 2);
                end

                % Create a new figure for each transition
                figure;



                % Create subplot for heatmap
                subplot(5, 1, 1:4)

                % Plot heatmap
                imagesc(transition_matrix);

                colorbar;

                % Create title string with state names
                transition_numbers = strsplit(char(transition_type), ' to ');

                % Change 'cwake' back to numeric representation for indexing scoreToState
                transition_numbers = cellfun(@(x) strrep(x, 'cwake', '1'), transition_numbers, 'UniformOutput', false);

                % Convert to numeric
                transition_numbers = cellfun(@str2num, transition_numbers);

                transition_states = scoreToState(transition_numbers);
                title_string = ['Heatmap for Transition Type: ' strjoin(transition_states, ' to ')];

                box off;

                xticks(1:fs_fiber*5:size(transition_matrix,2));
                xticklabels(-10:5:10);

                % Set title
                title(title_string);

                % Set x-axis label
                xlabel('Time (s)');

                % Set y-axis label
                ylabel('Transitions');



                % Groupby mouse
                transition_table = combined_transitions(transition_rows,:);

                % Unique mice
                unique_mouse = unique(transition_table.mouse);

                % Take average
                for m = 1:numel(unique_mouse)
                    indv_data = transition_table.Transition(transition_table.mouse == unique_mouse(m));
                    if size(indv_data,1) == 1
                        indv_matrix = indv_data{:};
                    else

                        indv_matrix = horzcat(indv_data{:});
                        indv_matrix = indv_matrix';
                    end

                    if size(indv_matrix,2)==1
                        m_indv_matrix = indv_matrix';
                    else
                        m_indv_matrix = mean(indv_matrix);
                    end

                    if m==1
                        mean_resps = m_indv_matrix;
                    else
                        mean_resps = [mean_resps ; m_indv_matrix];
                    end
                end


                % Create subplot for average line plot with SEM
                subplot(5, 1, 5)

                % Plot average line with SEM
                shadedErrorBar(1:length(mean_resps),mean(mean_resps),std(mean_resps)./sqrt(numel(unique_mouse)),'lineprops','r')
                hold on
                plot(mean_resps','k')
                hold off
                xlim([0,size(transition_matrix,2)]);
                xticks(0:fs_fiber*5:size(transition_matrix,2));
                xticklabels(-10:5:10);


                exportgraphics(gcf,fullfile(parent_folder,strcat(strjoin(transition_states, ' to '),'colorbar.pdf')),'Resolution',300,'ContentType','vector');


            end
        end


    end

end
