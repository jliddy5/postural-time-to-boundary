function [mean_ttb, med_ttb, min_ttb] = ttbBoundary(ttb_bound, n_min, plot_flag, labels)
%TTBBOUNDARY Calculates mean, median, and minimum TtB for each boundary.
%
% ARGUMENTS
% ttb_bound - Time-to-Boundary matrix for each boundary (n_samples x n_boundaries).
% Each column contains the TtB time series for one boundary. This function
% is designed to take the ttb_bound variable returned by the ttb function.
% 
% n_min - Number of minima to include in the estimation of min_ttb (positive integer).
% The default is 10% of the time series length.
%
% plot_flag - Boolean variable to request plots (0 = no plots, 1 = plots). Default is 0.
%
% labels - Cell array of boundary labels for making plots (required if plot_flag = 1).
%
% RETURNS
% mean_ttb - Vector containing the mean TtB for each boundary (excluding
% NaN values). Returns n x 1 vector with one value per boundary.
%
% med_ttb - Vector containing the median TtB for each boundary (excluding
% zeros). Returns n x 1 vector with one value per boundary.
%
% min_ttb - Vector containing the minimum TtB for each boundary (excluding
% zeros). Returns n x 1 vector with one value per boundary. The minimum is
% determined by taking the mean of a user specified number (n_min) of
% minima.
%
% ========================================================================%

%% Validation
arguments
    ttb_bound (:,:) double {mustBeNumeric, mustBeNonempty}
    n_min (1,1) double {mustBePositive, mustBeInteger} = NaN
    plot_flag (1,1) {mustBeNumericOrLogical} = 0
    labels (:,1) cell = {}
end

% Get dimensions
[n_samples, n_boundaries] = size(ttb_bound);

% Set default n_min if not provided
if isnan(n_min)
    n_min = round(0.1 * n_samples);
end

% Validate n_min is not larger than available data
if n_min > n_samples
    error('n_min (%d) cannot exceed the number of samples (%d).', n_min, n_samples);
end

% Validate labels if plotting
if plot_flag && isempty(labels)
    error('labels must be provided when plot_flag = 1.');
end

if plot_flag && length(labels) ~= n_boundaries
    error('Number of labels (%d) must match number of boundaries (%d).', length(labels), n_boundaries);
end

%% Compute statistics

% Preallocate output vectors
mean_ttb = zeros(n_boundaries, 1);
med_ttb = zeros(n_boundaries, 1);
min_ttb = zeros(n_boundaries, 1);

for i = 1:n_boundaries
    mean_ttb(i, 1) = mean(ttb_bound(:, i), 'omitnan');
    med_ttb(i, 1) = median(ttb_bound(ttb_bound(:, i) > 0, i), 'omitnan');
    min_ttb(i, 1) = ttbMinimum(ttb_bound(:, i), n_min);
end

%% Generate plots

if plot_flag
    
    % Mean TtB
    f = figure();
    theme(f, 'light');
    title('Mean TtB');
    bar(mean_ttb, 'FaceColor', [.5 .5 .5], 'EdgeColor', [0 0 0], 'LineWidth', 1.0);
    set(gca, 'XTick', [1:n_boundaries]', 'XTickLabel', labels);
    xlabel('Boundary');
    ylabel('Mean TtB (s)');
    
    % Median TtB
    f = figure();
    theme(f, 'light');
    title('Median TtB');
    bar(med_ttb, 'FaceColor', [.5 .5 .5], 'EdgeColor', [0 0 0], 'LineWidth', 1.0);
    set(gca, 'XTick', [1:n_boundaries]', 'XTickLabel', labels);
    xlabel('Boundary');
    ylabel('Median TtB (s)');
    
    % Minimum TtB
    f = figure();
    theme(f, 'light');
    title('Minimum TtB');
    bar(min_ttb, 'FaceColor', [.5 .5 .5], 'EdgeColor', [0 0 0], 'LineWidth', 1.0);
    set(gca, 'XTick', [1:n_boundaries]', 'XTickLabel', labels);
    xlabel('Boundary');
    ylabel('Minimum TtB (s)');
    
end

end
