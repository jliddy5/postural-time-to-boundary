%Clear workspace and command window.
clear; clc;

% Load example data
load('example_cop.mat');

%CONTENTS OF .MAT FILE ===================================================%
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
%   cop_y - center-of-pressure coordinate in the Ap direction
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
figure('Color', 'white');
subplot(2,1,1);
plot(time, cop_x, 'k');
xlabel('Time (s)');
ylabel('ML Position (mm)');
subplot(2,1,2);
plot(time, cop_y, 'k');
xlabel('Time (s)');
ylabel('AP Position (mm)');

% Force plate coordinates
x_fp = [fp_pts(:,1); fp_pts(1,1)];
y_fp = [fp_pts(:,2); fp_pts(1,2)];

% Base of support coordinates
x_bos = [bound_pts(:,1); bound_pts(1,1)];
y_bos = [bound_pts(:,2); bound_pts(1,2)];

% Plot of CoP within the base of support
figure('Color', 'white');
plot(cop_x, cop_y, 'k');
hold on;
plot(x_fp, y_fp, 'k');
scatter(fp_pts(:,1), fp_pts(:,2), 20, 'k', 'filled');
plot(x_bos, y_bos, 'k');
scatter(bound_pts(:,1), bound_pts(:,2), 20, 'k', 'filled');
hold off;
xlabel('ML Position (mm)');
ylabel('AP Position (mm)');
axis([(fp_pts(3,1) - 32) (fp_pts(1,1) + 32) (fp_pts(3,2) - 10) (fp_pts(1,2) + 10)]);

%Calculate Time-to-Contact ===============================================%

% Boundary Labels
bound_labels = {'F','FR','BR','B','BL','FL'};

% Sampling rate (Hz)
fs = 120;

% TtC method (1 = Riccio, 2 = Slobounov, 3 = Jerk)
extrapMethod = 2;

% Compute TtC
[ttc, ttc_bound, bound_crossed, bound_percent] = timetocontact(cop_x, cop_y, 1/fs, bound_pts, extrapMethod);

% Plot subset of virtual trajectories to check that method is working
plotVirtualTrajectory(cop_x, cop_y, 1/fs, x_bos, y_bos, ttc, extrapMethod);

% Plot TtC time series
figure('color','white');
plot(time, ttc, 'k');
xlabel('Time (s)');
ylabel('Time-To-Contact (s)');

% Plot TtC histogram with log normal fit
fitparm = lognfit(ttc);
logpdf = lognpdf([0:.01:max(ttc)]',fitparm(1),fitparm(2));
figure('color','white');
histogram(ttc,'Normalization','pdf','FaceColor',[50 50 50]./255,'EdgeColor',[0 0 0]);
hold on;
plot([0:.01:max(ttc)]',logpdf,'r');
xlabel('Time-to-Contact (s)');
ylabel('Probability density function');

% Compute mean, median, and minimum TtC for each boundary.
nMin = 50;
[ mean_ttc, med_ttc, min_ttc ] = ttcBoundary(ttc_bound, nMin, extrapMethod, bound_labels);