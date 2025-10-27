function [mean_ttb, med_ttb, min_ttb] = ttbBoundary(ttb, n_min, plot_flag, labels)
%TTBBOUNDARY Calculates Mean, Median, and Minimum TtB for each boundary.
%
% ARGUMENTS
% ttb - Time-to-Boundary for each boundary. This argument can either be an m
% x 1 vector or an m x n matrix with TtB values for n boundaries where m is
% the number of data points in the time series. This function was designed
% to take the ttb_bound variable returned by the timeToBoundary function.
% 
% n_min - Number of minima to include in the estimation of min_ttb. The
% default is 10% of the time series length.
%
% plot_flag - Boolean variable to request plots (0 = no plots, 1 = plots)
%
% labels - List of boundary labels for making plots.
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

% Default plot setting is off
if nargin == 2
    plot_flag = 0;
end

% Length of the time series and number of boundaries
[n_samples, n_boundaries] = size(ttb);

% Default number of minima is 10% of the time series length
if nargin == 1
    n_min = round(0.1 * n_samples);
    plot_flag = 0;
end

% Create mean, median, and minimum TtB vectors
mean_ttb = zeros(n_boundaries, 1);
med_ttb = zeros(n_boundaries, 1);
min_ttb = zeros(n_boundaries, 1);

for i = 1:n_boundaries
    mean_ttb(i, 1) = mean(ttb(:, i), 'omitnan');
    med_ttb(i, 1) = median(ttb(ttb(:, i) > 0, i), 'omitnan');
    min_ttb(i, 1) = ttbMinimum(ttb(:, i), n_min);
end

% Make Plots =============================================================%

if plot_flag
    
    % Mean TtB Bar Graph
    f = figure();
    theme(f, 'light');
    title('Mean TtB');
    bar(mean_ttb, 'FaceColor', [.5 .5 .5], 'EdgeColor', [0 0 0], 'LineWidth', 1.0);
    set(gca, 'XTick', [1:n_boundaries]', 'XTickLabel', labels);
    xlabel('Boundary');
    ylabel('Mean TtB (s)');
    
    % Median TtB Bar Graph
    f = figure();
    theme(f, 'light');
    title('Median TtB');
    bar(med_ttb, 'FaceColor', [.5 .5 .5], 'EdgeColor', [0 0 0], 'LineWidth', 1.0);
    set(gca, 'XTick', [1:n_boundaries]', 'XTickLabel', labels);
    xlabel('Boundary');
    ylabel('Median TtB (s)');
    
    % Minimum TtB Bar Graph
    f = figure();
    theme(f, 'light');
    title('Minimum TtB');
    bar(min_ttb, 'FaceColor', [.5 .5 .5], 'EdgeColor', [0 0 0], 'LineWidth', 1.0);
    set(gca, 'XTick', [1:n_boundaries]', 'XTickLabel', labels);
    xlabel('Boundary');
    ylabel('Minimum TtB (s)');
    
end

end
