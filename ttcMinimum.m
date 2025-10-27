function [ min_ttc ] = ttcMinimum( ttc, nMin )
%Estimates the minimum TtC as the average of the nMin lowest values in
% the TtC time series.

% Sort and compute the mean of the first nMin values.
sort_ttc = sort(ttc);
min_ttc = mean(sort_ttc(1:nMin),'omitnan');

end