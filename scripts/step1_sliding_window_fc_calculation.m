%% 00 Generate sliding window functional connectivity for all subject scans

%produces per subject per scan a region x region x window cell structure
%So it needs a cell of sub x ses (timeseries x roi)
%the idea for this workflow is that it could be used for a variable number
%of roi x roi per scan. For example, ICA-defined networks with a variable
%number of constituent regions (so not a 200 x 200 standard parcellation

%% Init
close all
clear all

%% Path
%main
paths.home = (cd);
addpath(genpath(paths.home))

paths.out = fullfile(paths.home,'results','dynamic_flex');
if ~exist(paths.out)
    mkdir(paths.out);
end

%% Data
load(fullfile(paths.home,'data','ROI_timeseries_PSC_noGSR_scrub.mat'));
%a cell structure of sub x ses x timepoint x roi

%example
time_series_data = roi_timeseries.schaefer232;

%% Define parameters
tr_seconds = 1.4; %for note keeping
window_length = 30; % Number of timepoints per window (42s)
overlap_length = window_length / 2; % 50% overlap
total_timepoints = 516;
num_windows = floor((total_timepoints - window_length) / overlap_length) + 1; %total number of windows
n_sub = size(time_series_data,1);
n_ses = size(time_series_data,2);

%% Begin
% Initialize an empty cell array to store functional connectivity matrices for each window
dynamic_correlation_matrices = cell(n_sub,n_ses);

% Iterate over scans
for sub = 1:n_sub
    for ses = 1:n_ses
        scan_data = time_series_data{sub,ses}'; % Time series data for the current scan
        num_rois = size(scan_data, 1); % Number of ROIs for this scan

        % Initialize an empty array to store functional connectivity matrices for this scan
        scan_functional_matrices = zeros(num_rois, num_rois, num_windows);

        % Iterate over windows
        for window_start = 1:overlap_length:(total_timepoints - window_length + 1)
            window_end = window_start + window_length - 1;

            % Extract data for the current window
            window_data = scan_data(:, window_start:window_end);

            % Calculate the functional connectivity matrix for the window
            correlation_matrix = corrcoef(window_data');

            % Store the functional connectivity matrix
            window_index = ceil(window_start / overlap_length);
            scan_functional_matrices(:,:,window_index) = correlation_matrix;
        end

        % Store the functional connectivity matrices for this scan in the cell array
        dynamic_correlation_matrices{sub,ses} = scan_functional_matrices;
    end
end
% Now, dynamic_correlation_matrices is a cell array containing functional
% connectivity matrices for each scan with a variable number of region per
% window

%% Save output

save(fullfile(paths.out,'sliding_window_fc.mat'),"dynamic_correlation_matrices")
