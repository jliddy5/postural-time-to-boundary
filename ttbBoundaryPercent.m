function [boundary_percent] = ttbBoundaryPercent(which_boundary, n_boundaries)
%TTBBOUNDARYPERCENT Calculates the percentage of virtual contacts attributed to each boundary.
% 
% ARGUMENTS
% which_boundary - Vector of crossed boundaries for each instant in time (integers 1 to n_boundaries)
% n_boundaries - Total number of boundaries (positive integer)
%
% RETURNS
% boundary_percent - Vector of size n_boundaries x 1 containing the 
% percentage of contacts for each boundary.
%
%=========================================================================%

%% Validation
arguments
    which_boundary (:,1) double {mustBeNumeric, mustBeNonempty, mustBeInteger}
    n_boundaries (1,1) double {mustBePositive, mustBeInteger}
end

% Validate boundary indices are in valid range
if any(which_boundary < 1) || any(which_boundary > n_boundaries)
    error('which_boundary values must be integers between 1 and %d.', n_boundaries);
end

%% Compute percentages

% Number of data points
n_samples = length(which_boundary);

% Preallocate percentage vector
boundary_percent = zeros(n_boundaries, 1);

% Compute percentage for each boundary
for i = 1:n_boundaries
    boundary_percent(i, 1) = (sum(which_boundary == i) / n_samples) * 100;
end

end
