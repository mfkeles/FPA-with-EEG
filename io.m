classdef io
    properties
        root_dir % Root directory to be searched
    end

    methods
        % Constructor
        function obj = io(root_dir)
            obj.root_dir = root_dir; % Store the root directory
        end

        % Function to find all directories starting with "animal"
        function mouse_dirs = find_mouse_dirs(obj)
            all_dirs = dir(fullfile(obj.root_dir, 'animal*')); % Search for directories starting with "animal"
            mouse_dirs = all_dirs([all_dirs.isdir]); % Filter out non-directory entries
        end

        % Function to find all trial directories in a mouse directory
        function trial_dirs = find_trial_dirs(obj, mouse_dir)
            pattern = fullfile(obj.root_dir, mouse_dir.name, '*ZT*-*'); % Define search pattern
            trial_dirs = dir(pattern); % Search for trial directories
        end
    end

    methods (Static)
        % Static function to find the required CSV, EDF, and TSV files in a trial directory
        function files = find_files(trial_dir)
            % Find all CSV files
            csv_files = dir(fullfile(trial_dir.folder, trial_dir.name, '*.csv'));

            % Find the smallest CSV file in terms of size
            if ~isempty(csv_files)
                [~, min_idx] = min([csv_files.bytes]);
                csv_file = csv_files(min_idx);
            else
                csv_file = [];
            end

            % Find the required EDF and TSV files
            edf_file = dir(fullfile(trial_dir.folder, trial_dir.name, '*export.edf'));
            tsv_file = dir(fullfile(trial_dir.folder, trial_dir.name, '*.tsv'));

            % Store the files in a structure
            files = struct('csv_file', csv_file, 'edf_file', edf_file, 'tsv_file', tsv_file);
        end
    end

    methods
        % Function to create the output table
        function table_out = create_table(obj)
            mouse_dirs = obj.find_mouse_dirs(); % Find all mouse directories

            table_data = cell(0, 5); % Initialize an empty cell array for table data
            for i = 1:length(mouse_dirs)
                trial_dirs = obj.find_trial_dirs(mouse_dirs(i)); % Find all trial directories for the current mouse directory
                for j = 1:length(trial_dirs)
                    files = io.find_files(trial_dirs(j)); % Find the required files for the current trial directory
                    % Check if all required files are present
                    if ~isempty(files.csv_file) && ~isempty(files.edf_file) && ~isempty(files.tsv_file)
                        % Add the data to the table_data cell array
                        table_data = [table_data; {mouse_dirs(i).name, trial_dirs(j).name, ...
                            fullfile(files.csv_file.folder, files.csv_file.name), ...
                            fullfile(files.edf_file.folder, files.edf_file.name), ...
                            fullfile(files.tsv_file.folder, files.tsv_file.name)}];
                    end
                end
            end

            % Convert the cell array to a table with appropriate column names
            table_out = array2table(table_data, 'VariableNames', {'mouse_dir', 'trial_dir', 'fiber_dir', 'eeg_dir', 'score_dir'});
        end
    end

    methods (Static)
        function save_processed_data(fiber_dir, fiber_chunks_with_scores, eeg_chunks_table)
            % Get the parent folder of the fiber file
            parent_folder = fileparts(fiber_dir);

            % Create a new subfolder for the processed data if it doesn't exist
            processed_data_folder = fullfile(parent_folder, 'processed_data');
            if ~exist(processed_data_folder, 'dir')
                mkdir(processed_data_folder);
            end

            % Save fiber_chunks_with_scores
            fiber_chunks_file = fullfile(processed_data_folder, 'fiber_chunks_with_scores.mat');
            save(fiber_chunks_file, 'fiber_chunks_with_scores');

            % Save eeg_chunks_table
            eeg_chunks_file = fullfile(processed_data_folder, 'eeg_chunks_table.mat');
            save(eeg_chunks_file, 'eeg_chunks_table');
        end
        function zt_portion = extract_zt(row)
            % Extract the cell content from the "trial_dir" column
            trial_dir = row.trial_dir{1};

            % Find the position of the last dash
            last_dash_pos = find(trial_dir == '-', 1, 'last');

            % Extract the portion starting with 'ZT'
            zt_portion = trial_dir(last_dash_pos + 1:end);
        end
    end
end
