function [ mean_ttc, med_ttc, min_ttc ] = ttcBoundary( ttc, nMin, PLOT, labels )
%TTC_BOUND Calculates Mean, Medan, and Minimum TtC for each boundary.
%
% ARGUMENTS
% ttc - Time-to-Contact for each boundary. This argument can either be an m
% x 1 vector or an m x n matrix with TtC values for n boundaries where m is
% the number of data points in the time series. This function was designed
% to take the bound_ttc variable returned by the timetocontact function.
% 
% nMin - Number of minima to include in the estimation of min_ttc. The
% default is 10% of the time series length.
%
% PLOT - boolean variable to request plots (0 = no plots, 1 = plots)
%
% labels - List of boundary labels for making plots.
%
% RETURNS
% mean_ttc - Vector containing the mean TtC for each boundary (excluding
% zeros). Returns n x 1 vector with one value per boundary.
%
% med_ttc - Vector containing the median TtC for each boundary (excluding
% zeros). Returns n x 1 vector with one value per boundary.
%
% min_ttc - Vector containing the minimum TtC for each boundary (excluding
% zeros). Returns n x 1 vector with one value per boundary. The minimum is
% determined by taking the mean of a user specified number (n_mins) of
% minima.
%
% ========================================================================%

% Default plot setting is off.
if nargin == 2
    PLOT = 0;
end

% Length of the time series and number of boundaries.
[N, nBoundaries] = size(ttc);

%default number of minima is 10% of the time series length
if nargin == 1
    nMin = round(0.1 * N);
    PLOT = 0;
end

% Create mean, median, and minimum TtC vectors
mean_ttc = zeros(nBoundaries,1);
med_ttc = zeros(nBoundaries,1);
min_ttc = zeros(nBoundaries,1);

for i = 1:nBoundaries
    mean_ttc(i,1) = mean(ttc(:,i),'omitnan');
    med_ttc(i,1) = median(ttc(ttc(:,i) > 0,i),'omitnan');
    min_ttc(i,1) = ttcMinimum(ttc(:,i),nMin);
end

%Make Plots ==============================================================%

if PLOT
    
    %Mean TtC Bar Graph
    figure('Color', 'white');
    title('Mean TtC');
    bar(mean_ttc,'FaceColor',[.5 .5 .5],'EdgeColor',[0 0 0], 'LineWidth', 1.0);
    set(gca,'XTick', [1:nBoundaries]', 'XTickLabel', labels);
    xlabel('Boundary');
    ylabel('Mean TtC(s)');
    
    %Median TtC Bar Graph
    figure('Color', 'white');
    title('Median TtC');
    bar(med_ttc,'FaceColor',[.5 .5 .5],'EdgeColor',[0 0 0], 'LineWidth', 1.0);
    set(gca,'XTick', [1:nBoundaries]', 'XTickLabel', labels);
    xlabel('Boundary');
    ylabel('Median TtC(s)');
    
    %Minimum TtC Bar Graph
    figure('Color', 'white');
    title('Minimum TtC');
    bar(min_ttc,'FaceColor',[.5 .5 .5],'EdgeColor',[0 0 0], 'LineWidth', 1.0);
    set(gca,'XTick', [1:nBoundaries]', 'XTickLabel', labels);
    xlabel('Boundary');
    ylabel('Minimum TtC(s)');
    
end

end