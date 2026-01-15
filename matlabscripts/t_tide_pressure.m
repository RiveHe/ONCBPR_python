% Load the data
addpath('F:\aloha\ALOHA Cabled Observatory Database\t_tide_v1.5beta');
savepath;

filePath = 'prs_pressure_2014_raw_pressure_only.csv';

opts = detectImportOptions(filePath, 'Delimiter', ',', 'NumHeaderLines', 0);
data = readtable(filePath, opts);

% Convert Time_UTC to datetime format if not already
if ~isdatetime(data.Time_UTC)
    data.Time_UTC = datetime(data.Time_UTC, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC');
end

% Calculate the difference between consecutive timestamps
time_diffs = diff(data.Time_UTC);

% Define a gap threshold, e.g., more than 1 hour
gap_threshold = hours(1);

% Find indices where the gap exceeds the threshold
gap_indices = find(time_diffs > gap_threshold);

% Start and end indices of continuous segments
segment_starts = [1; gap_indices + 1];
segment_ends = [gap_indices; height(data)];

% Define a limit for sub-segmentation based on the hourly data
year_limit_hours = 8760;  % Normal year - for leap year adjust to 8784

% Pre-process segments to fit within the year-long limit in hours
final_segments = [];
for i = 1:length(segment_starts)
    current_start = segment_starts(i);
    while current_start <= segment_ends(i)
        % Calculate the end of the current sub-segment within year limit in hours
        current_end = min(current_start + year_limit_hours - 1, segment_ends(i));
        final_segments = [final_segments; current_start, current_end];
        current_start = current_end + 1;  % Move to the next segment start
    end
end

% Initialize results table
results = table;

% Tidal analysis and prediction for each pre-segmented interval
for i = 1:size(final_segments, 1)
    segment = data(final_segments(i, 1):final_segments(i, 2), :);
    t_segment = datenum(segment.Time_UTC);
    mean_pressure = mean(segment.RawPressure);
    pressure_anomaly = segment.RawPressure - mean_pressure;

    % Tidal analysis using t_tide
    [nameu, fu, tidecon, xout] = t_tide(pressure_anomaly, 'start', ...
        t_segment(1), 'latitude', 47.958352, 'error', 'cboot', 'lsq', 'normal');

    predicted_tides = t_predic(t_segment, nameu, fu, tidecon, 'latitude', 47.948583);

    detided_data = pressure_anomaly - predicted_tides;

    % Second round of t_tide analysis
    [nameu, fu, tidecon, xout] = t_tide(detided_data, 'start', ...
        t_segment(1), 'latitude', 47.958352, 'error', 'cboot', 'lsq', 'normal');

    predicted_tides = t_predic(t_segment, nameu, fu, tidecon, 'latitude', 47.948583);

    detided_data2 = detided_data - predicted_tides;

    % Plotting the results
    figure;
    hold on;
    plot(segment.Time_UTC, pressure_anomaly, 'b', 'DisplayName', 'Original Data');
    plot(segment.Time_UTC, detided_data2, 'g', 'DisplayName', 'Detided Data');
    plot(segment.Time_UTC, predicted_tides, 'r', 'DisplayName', 'Predicted Tides');

    xlabel('Time');
    ylabel('Pressure (decibar)');
    title(sprintf('Tidal Analysis and Detiding for Segment %d', i));
    legend show;
    grid on;

    % Save results for this segment
    segment_results = table(segment.Time_UTC, segment.RawPressure, detided_data2, ...
        'VariableNames', {'Time_UTC', 'Raw Pressure', 't_tide_detided'});
    results = [results; segment_results]; 
end

% Save the final results to a CSV file
writetable(results, 'prs_pressure_detided.csv');
%%
%% load the detided data

addpath('F:\aloha\ALOHA Cabled Observatory Database');
savepath;

% Load the detided data
filePathDetided = 'prs_pressure_detided.csv';
optsDetided = detectImportOptions(filePathDetided, 'Delimiter', ',', 'NumHeaderLines', 0, 'VariableNamingRule', 'preserve');
detidedData = readtable(filePathDetided, optsDetided);

% Convert Time_UTC to datetime format if not already
if ~isdatetime(detidedData.Time_UTC)
    detidedData.Time_UTC = datetime(detidedData.Time_UTC, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC');
end

%% Check for and remove any NaN or Inf values from the detided pressure data
validIndices = isfinite(detidedData.t_tide_detided);
filteredTime_UTC = detidedData.Time_UTC(validIndices);
filteredPressureData = detidedData.t_tide_detided(validIndices);

% Define the sampling frequency
fs = 2; % Sampling rate of 2 samples per second

% Convert the desired cutoff period (3 days) into a frequency
cutoff_period = 3 * 24 * 3600; % 3 days in seconds
cutoff_frequency = 1 / cutoff_period; % Cutoff frequency in Hz

% Design a low-pass filter for the cutoff frequency
d = designfilt('lowpassiir', 'FilterOrder', 8, ...
               'HalfPowerFrequency', cutoff_frequency, ...
               'SampleRate', fs);

% Apply the filter to the detided pressure data
filteredPressure = filtfilt(d, filteredPressureData);

% Create a table with the filtered data
filteredTable = table(filteredTime_UTC, filteredPressure, ...
    'VariableNames', {'Time_UTC', 'Filtered_Pressure'});

% Save the filtered data to a CSV file
writetable(filteredTable, 'prs_pressure_detided_filtered.csv');

%% Plot the results
figure;
%plot(filteredTime_UTC, filteredPressureData, 'g', 'DisplayName', 'Detided Pressure');
hold on;

% Plot the filtered pressure data
plot(filteredTime_UTC, filteredPressure, 'r', 'DisplayName', 'Filtered Detided Pressure (Low-pass)');

xlabel('Time');
ylabel('Pressure (decibar)');
title('Low-pass Filtered Detided Pressure Data (Periods > 3 days)');
legend show;
grid on;
hold off;



%%
function outliers = rolling_mad_outlier_detection(series, window_size, thresh)
    if nargin < 3
        thresh = 2.6; % Default threshold if not provided
    end
    
    % Ensure the series is a column vector
    series = series(:);
    
    % Calculate the rolling median
    rolling_median = movmedian(series, window_size, 'Endpoints', 'shrink');
    
    % Calculate the rolling MAD (Median Absolute Deviation)
    rolling_mad = movmad(series, window_size, 'Endpoints', 'shrink');
    
    % Calculate the modified Z-score
    median_absolute_deviation = 1.4826;
    modified_z_score = abs(series - rolling_median) ./ (median_absolute_deviation * rolling_mad);
    
    % Identify outliers
    outliers = modified_z_score > thresh;
end

% Example usage with filteredPressure data:
% Define the series and parameters
series = filteredPressure; % Replace with your data
window_size = 200; % Define the window size (adjust as needed)

% Detect outliers
outliers = rolling_mad_outlier_detection(series, window_size);

% Remove outliers from the series
cleaned_series = series;
cleaned_series(outliers) = NaN; % Set outliers to NaN

% Optionally, you can perform further smoothing on the cleaned series
rolling_window = 200; % Adjust the rolling window size as needed
smoothed_series = movmean(cleaned_series, rolling_window, 'omitnan');

% Plotting the original, cleaned, and smoothed series
figure;
plot(series, 'k', 'DisplayName', 'Original Series');
hold on;
plot(cleaned_series, 'b', 'DisplayName', 'Series without Outliers');
plot(smoothed_series, 'r', 'DisplayName', 'Smoothed Series');
xlabel('Time');
ylabel('Pressure (decibar)');
title('Outlier Removal and Smoothing');
legend show;
grid on;
hold off;
%%

