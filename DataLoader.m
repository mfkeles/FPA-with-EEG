classdef DataLoader
    properties
        data_table % Input data table
        current_row % Current row number

        score_data % Loaded data from the score_dir file
        eeg_data % Loaded data from the eeg_dir file
        fiber_data % Loaded data from the fiber_dir file
        
    end

    methods
        % Constructor
        function obj = DataLoader(data_table)
            obj.data_table = data_table;
            obj.current_row = 1;
        end

        % Function to load a specific row
        function obj = load_row(obj, row_number)
            if row_number > 0 && row_number <= height(obj.data_table)
                obj.current_row = row_number;
                obj = obj.load_score_data();
                obj = obj.load_eeg_data();
                obj = obj.load_fiber_data();
            else
                error('Invalid row number.');
            end
        end

        function obj = load_score_data(obj)
            score_file = obj.data_table.score_dir{obj.current_row};


            % Read the first 20 lines of the file
            fid = fopen(score_file);
            lines = cell(20, 1);
            for i = 1:20
                lines{i} = fgetl(fid);
            end
            fclose(fid);

            % Find the line starting with "Date"
            date_line_idx = find(cellfun(@(x) strncmp(x, 'Date', 4), lines), 1);
            if isempty(date_line_idx)
                error('Could not find a line starting with "Date".');
            end

            % Set NumHeaderLines to the line that starts with "Date"
            num_header_lines = date_line_idx - 1;

            % Load the score data with the correct NumHeaderLines
            temp_score_data = readtable(score_file, 'FileType', 'text', 'NumHeaderLines', num_header_lines);

            % Check if there are 6 columns
            if width(temp_score_data) ~= 6
                error('The score data must have exactly 6 columns.');
            end

            % Remove the 6th column
            temp_score_data(:, 6) = [];

            % Check and correct unique values
            unique_values = unique(temp_score_data{:, end});
            invalid_values = setdiff(unique_values, [1, 2, 3, 4, 255]);
            if ~isempty(invalid_values)
                for i = 1:numel(invalid_values)
                    temp_score_data{temp_score_data{:, end} == invalid_values(i), end} = 255;
                end
            end

            obj.score_data = temp_score_data;
        end
        % Function to load the EEG data from the eeg_dir file
        function obj = load_eeg_data(obj)
            eeg_file = obj.data_table.eeg_dir{obj.current_row};
            obj.eeg_data = edfread(eeg_file);
        end

        % Function to load the fiber data from the fiber_dir file
        function obj = load_fiber_data(obj)
            fiber_file = obj.data_table.fiber_dir{obj.current_row};
            obj.fiber_data = readtable(fiber_file);
        end
    end
end
