% Test script for time-to-boundary calculations
% Clear workspace and command window
clear; clc;

% Load example data
load('example_cop.mat');

% CONTENTS OF .MAT FILE ===================================================%
% 30 s of quiet standing collected @ 120 Hz with AMTI OR6-7 Force Platform
% providing the three components of the ground reaction force (Fx, Fy, &
% Fz) and net moment (Mx, My, & Mz).
%
% Center-of-pressure was calculated after filtering the raw forces and
% moments using an 8th order Butterworth, low-pass filter with a cutoff
% frequency of 20 Hz. The force plate dimensions are 464 mm (width) x 508
% mm (length). The force plate origin is located at the bottom right hand
% corner with positive x and y axes extending to the right and forward,
% respectively. Thus, all CoP coordinates are positive and bounded on the
% intervals [0,464] and [0,508] in the ML and AP directions, respectively. 
% 
% For this example the anatomical landmarks selected to define the
% base-of-support were the 2nd phalanx of the Hallux, 5th metatarsal, and
% the most dorsal point of the calcaneus for each foot.
%
% Variables
%   cop_x - center-of-pressure coordinate in the ML direction
%   cop_y - center-of-pressure coordinate in the AP direction
%   time - time series with dt = 1/samp_rate = 1/120
%   bound_pts - coordinates for the landmarks selected to define the 
%               boundaries of the base-of-support (column 1 = x coordinate,
%               column 2 = y coordinate). The landmarks are in the following
%               order: right toe, right m5, right heel, left heel, left m5,
%               and left toe. ** BOUNDARY LANDMARKS MUST BE ORGANIZED IN
%               CLOCKWISE FASHION, BUT MAY START WITH ANY LANDMARK.
%   fp_pts - coordinates of the corners of the force platform starting with
%            at the top right and moving clockwise (column 1 = x coordinate,
%            column 2 = y coordinate).
% ========================================================================%

% Create CoP plots =======================================================%
% Plot of CoP components
f = figure();
theme(f, 'light');
subplot(2, 1, 1);
plot(time, cop_x, 'k');
xlabel('Time (s)');
ylabel('ML Position (mm)');
subplot(2, 1, 2);
plot(time, cop_y, 'k');
xlabel('Time (s)');
ylabel('AP Position (mm)');

% Force plate coordinates
x_fp = [fp_pts(:, 1); fp_pts(1, 1)];
y_fp = [fp_pts(:, 2); fp_pts(1, 2)];

% Base of support coordinates
x_bos = [bound_pts(:, 1); bound_pts(1, 1)];
y_bos = [bound_pts(:, 2); bound_pts(1, 2)];

% Plot of CoP within the base of support
f = figure();
theme(f, 'light');
plot(cop_x, cop_y, 'k');
hold on;
plot(x_fp, y_fp, 'k');
scatter(fp_pts(:, 1), fp_pts(:, 2), 20, 'k', 'filled');
plot(x_bos, y_bos, 'k');
scatter(bound_pts(:, 1), bound_pts(:, 2), 20, 'k', 'filled');
hold off;
xlabel('ML Position (mm)');
ylabel('AP Position (mm)');
axis([(fp_pts(3, 1) - 32) (fp_pts(1, 1) + 32) (fp_pts(3, 2) - 10) (fp_pts(1, 2) + 10)]);

% Calculate Time-to-Boundary =============================================%

% Boundary Labels
bound_labels = {'F', 'FR', 'BR', 'B', 'BL', 'FL'};

% Sampling rate (Hz)
fs = 120;

% TtB method (1 = Riccio, 2 = Slobounov, 3 = Jerk)
extrap_method = 2;

% Compute TtB
[ttb, ttb_bound, bound_crossed, bound_percent] = ttb(cop_x, cop_y, 1/fs, bound_pts, extrap_method);

% Plot TtB time series
f = figure();
theme(f, 'light');
plot(time, ttb, 'k');
xlabel('Time (s)');
ylabel('Time-to-Boundary (s)');

% Plot TtB histogram with log normal fit
fit_params = lognfit(ttb);
log_pdf = lognpdf([0:0.01:max(ttb)]', fit_params(1), fit_params(2));
f = figure();
theme(f, 'light');
histogram(ttb, 'Normalization', 'pdf', 'FaceColor', [50 50 50]./255, 'EdgeColor', [0 0 0]);
hold on;
plot([0:0.01:max(ttb)]', log_pdf, 'r', 'LineWidth', 2);
xlabel('Time-to-Boundary (s)');
ylabel('Probability Density Function');
legend('Data', 'Lognormal Fit', 'Location', 'best');

% Compute mean, median, and minimum TtB for each boundary
n_min = 50;
plot_flag = 1; % Set to 1 to display plots, 0 to suppress
[mean_ttb, med_ttb, min_ttb] = ttbBoundary(ttb_bound, n_min, plot_flag, bound_labels);

% Display summary statistics
fprintf('\n=== Time-to-Boundary Summary Statistics ===\n');
fprintf('Overall TtB: Mean = %.3f s, Median = %.3f s, Min = %.3f s\n', ...
    mean(ttb, 'omitnan'), median(ttb, 'omitnan'), min(ttb));
fprintf('\nBoundary-specific TtB:\n');
for i = 1:length(bound_labels)
    fprintf('  %s: Mean = %.3f s, Median = %.3f s, Min = %.3f s (%.1f%%)\n', ...
        bound_labels{i}, mean_ttb(i), med_ttb(i), min_ttb(i), bound_percent(i));
end
fprintf('\n');
