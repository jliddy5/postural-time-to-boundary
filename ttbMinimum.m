function [min_ttb] = ttbMinimum(ttb, n_min)
%TTBMINIMUM Estimates the minimum TtB as the average of the n_min lowest values
% in the TtB time series.
%
% ARGUMENTS
% ttb - Time-to-Boundary time series vector
% n_min - Number of minimum values to average
%
% RETURNS
% min_ttb - Average of the n_min lowest TtB values

% Sort and compute the mean of the first n_min values
sorted_ttb = sort(ttb);
min_ttb = mean(sorted_ttb(1:n_min), 'omitnan');

end
