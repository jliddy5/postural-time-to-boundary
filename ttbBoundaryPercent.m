function [boundary_percent] = ttbBoundaryPercent(which_boundary, n_boundaries)
%TTBBOUNDARYPERCENT Calculates the percentage of virtual contacts attributed to each boundary.
% 
% ARGUMENTS
% which_boundary - Vector of crossed boundaries for each instant in time
% n_boundaries - Number of boundaries
%
% RETURNS
% boundary_percent - Vector of size n_boundaries x 1 containing the 
% percentage of contacts for each boundary.
%
%=========================================================================%

% Number of data points
n_samples = length(which_boundary);

% Percentage vector
boundary_percent = zeros(n_boundaries, 1);

% Compute percentages
for i = 1:n_boundaries
    boundary_percent(i, 1) = (sum(which_boundary == i) / n_samples) * 100;
end

end
