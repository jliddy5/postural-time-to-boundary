function [ttb, ttb_bound, bound_crossed, bound_percent] = ttb(r_x, r_y, dt, bounds, extrap_method)
%TTB Calculates Time-to-Boundary (tau)
%
% ARGUMENTS
% r_x - ML position of the center of pressure or center of mass (numeric vector)
%
% r_y - AP position of the center of pressure or center of mass (numeric vector)
%
% dt - Unit change in time between samples (1/fs) (positive scalar)
%
% bounds - Matrix of boundary coordinates (x in the first column, y in the
% second). Boundary coordinates should be entered in ordered clockwise, but
% need not begin with any particular coordinate. If there are n boundary points (and
% therefore n boundaries) this matrix should be size n x 2 (minimum 3 points).
%
% extrap_method - Determines the method for estimating tau (1 or 2).
%          The default method is Slobounov (method 2)
%
%     Method 1 (Riccio)
%     Linear Equation:
%     A * tau + B = 0
%     where tau = virtual time-to-boundary
%       A = [v_y(t) - m_bound * v_x(t)]
%       B = [(p_y(t) - y_bound) - m_bound * (p_x(t) - x_bound)]
%     
%     Method 2 (Slobounov)
%     Quadratic Equation:
%     A * tau^2 + B * tau + C = 0
%     where tau = virtual time-to-boundary
%       A = [a_y(t) - m_bound * a_x(t)] / 2
%       B = [v_y(t) - m_bound * v_x(t)]
%       C = [(p_y(t) - y_bound) - m_bound * (p_x(t) - x_bound)]
%
% RETURNS
% ttb - Time series of minimum Time-to-Boundary in seconds. Minimum TtB is
% determined by taking smallest of the positive, real roots to the
% polynomial equations corresponding to the time it would take the
% extrapolated CoP/CoM trajectory to reach a boundary based on the
% instantaneous position, velocity, and acceleration.
%
% ttb_bound - Minimum Time-to-Boundary values for each boundary. The minimum
% of the positive, real roots for each point in time is selected. If no
% positive, real roots exist a value of NaN is returned.
%
% bound_crossed - Vector indicating which boundary had minimum TtB at each
% time point.
%
% bound_percent - Calculates the percentage of virtual contacts to each
% boundary. Uses the ttbBoundaryPercent function.
%
% ========================================================================%

%% Validation
arguments
    r_x (:,1) double {mustBeNumeric, mustBeNonempty}
    r_y (:,1) double {mustBeNumeric, mustBeNonempty}
    dt (1,1) double {mustBePositive}
    bounds (:,2) double {mustBeNumeric}
    extrap_method (1,1) double {mustBeMember(extrap_method, [1, 2])} = 2
end

% Additional validation: r_x and r_y must be same length
if length(r_x) ~= length(r_y)
    error('r_x and r_y must have the same length.');
end

% Additional validation: bounds must have at least 3 points
if size(bounds, 1) < 3
    error('bounds must have at least 3 boundary points to form a closed polygon.');
end

%% Compute derivatives

% Velocity
v_x = gradient(r_x, dt);
v_y = gradient(r_y, dt);

% Acceleration
a_x = gradient(v_x, dt);
a_y = gradient(v_y, dt);

%% Setup polynomial coefficients for each boundary

% Length of the data
n_samples = length(r_x);

% Number of boundaries
[n_boundaries, ~] = size(bounds);

% Extract boundary coordinates
x_bound = bounds(:, 1);
y_bound = bounds(:, 2);

% Compute boundary slopes (vectorized)
% Each boundary is defined as a line segment from the previous point to the current point.
% For the first boundary, connect from the last point (closing the polygon).
x_prev = [bounds(end, 1); bounds(1:end-1, 1)];
y_prev = [bounds(end, 2); bounds(1:end-1, 2)];

slope = (y_prev - y_bound) ./ (x_prev - x_bound);

% Handle vertical boundaries (infinite slope)
slope(isinf(slope)) = 1e16;

% Preallocate polynomial coefficient cells
coef_A = cell(n_boundaries, 1);
coef_B = cell(n_boundaries, 1);
coef_C = cell(n_boundaries, 1);

% Polynomial matrix
poly_coef = cell(n_boundaries, 1);

% Roots matrix
poly_roots = cell(n_boundaries, 1);

% For each boundary, compute polynomial coefficients
for i = 1:n_boundaries
    
    % Determine polynomial coefficients (A, B, C)
    coef_A{i} = (a_y - slope(i) .* a_x) ./ 2;
    coef_B{i} = (v_y - slope(i) .* v_x);
    coef_C{i} = (r_y - y_bound(i)) - slope(i) .* (r_x - x_bound(i));
    
    % Set A coefficient to 0 for linear method (Method 1)
    if extrap_method == 1
        coef_A{i} = zeros(n_samples, 1);
    end
    
    % Create polynomial and root matrices
    % For Method 1: [0 B C] → A·τ + B = 0 (linear)
    % For Method 2: [A B C] → A·τ² + B·τ + C = 0 (quadratic)
    poly_coef{i} = [coef_A{i} coef_B{i} coef_C{i}];
    poly_roots{i} = nan(n_samples, extrap_method);
    
end

%% Compute time-to-boundary

% Preallocate output arrays
ttb = zeros(n_samples, 1);
ttb_bound = zeros(n_samples, n_boundaries);
bound_crossed = zeros(n_samples, 1);
all_roots = zeros(n_samples, n_boundaries * extrap_method);

% For each time point
for t = 1:n_samples
      
    % For each boundary
    for i = 1:n_boundaries
        
        % Find roots of polynomial and add to full matrix
        poly_roots{i}(t, :) = roots(poly_coef{i}(t, :))';
        root_start_idx = extrap_method * (i - 1) + 1;
        root_end_idx = extrap_method * i;
        all_roots(t, root_start_idx:root_end_idx) = poly_roots{i}(t, :);
        
        % Find real roots greater than or equal to 0
        idx_real_pos = find(poly_roots{i}(t, :) >= 0 & imag(poly_roots{i}(t, :)) == 0);
        
        % Find minimum time-to-boundary for this boundary
        if isempty(idx_real_pos)
            ttb_bound(t, i) = NaN;
        else
            ttb_bound(t, i) = min(poly_roots{i}(t, idx_real_pos));
        end
    end
    
    % Find overall minimum TtB across all boundaries for this time point
    % Find real roots greater than or equal to 0
    idx_real_pos = find(all_roots(t, :) >= 0 & imag(all_roots(t, :)) == 0);
    
    % Find minimum time-to-boundary from the positive real roots
    [ttb(t), min_idx] = min(all_roots(t, idx_real_pos));
    
    % Identify the boundary that was crossed. Divides by the
    % number of roots per boundary and increments to the next 
    % integer to get the corresponding boundary number.
    bound_crossed(t) = ceil(idx_real_pos(min_idx) / extrap_method);
    
end

% Compute the percent of minimum time-to-boundary to each boundary
bound_percent = ttbBoundaryPercent(bound_crossed, n_boundaries);

end
