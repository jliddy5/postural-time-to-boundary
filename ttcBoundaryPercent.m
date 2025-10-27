function [ boundaryPercent ] = ttcBoundaryPercent( whichBoundary, nBoundaries )
% Calculates the percentage of virtual contacts attributed to each boundary.
% 
% ARGUMENTS
% counts - vector of crossed boundaries for each instant in time
% n_bounds - number of boundaries
%
% RETURNS
% boundaryPercent - Vector of size nBoundaries x 1 containing the 
% percentage of contacts for each boundary.
%
%=========================================================================%

% Number of data points
N = length(whichBoundary);

% Percentage vector
boundaryPercent = zeros(nBoundaries,1);

% Compute percentages
for i = 1:nBoundaries
    boundaryPercent(i,1) = (sum(whichBoundary == i) / N) * 100;
end

end