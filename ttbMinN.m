function [min_ttb] = ttbMinN(ttb, n_min)
%TTBMINN Estimates the minimum TtB as the average of the n_min lowest values
% in the TtB time series.
%
% ARGUMENTS
% ttb - Time-to-Boundary time series vector (numeric column vector)
% n_min - Number of minimum values to average (positive integer)
%
% RETURNS
% min_ttb - Average of the n_min lowest TtB values
%
%=========================================================================%

%% Validation
arguments
    ttb (:,1) double {mustBeNumeric, mustBeNonempty}
    n_min (1,1) double {mustBePositive, mustBeInteger}
end

% Validate n_min does not exceed vector length
if n_min > length(ttb)
    error('n_min (%d) cannot exceed the length of ttb vector (%d).', n_min, length(ttb));
end

%% Compute minimum

% Sort and compute the mean of the first n_min values
sorted_ttb = sort(ttb);
min_ttb = mean(sorted_ttb(1:n_min), 'omitnan');

end
