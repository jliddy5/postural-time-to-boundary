function [ttb, ttb_bound, bound_crossed, bound_percent] = timeToBoundary(r_x, r_y, dt, bounds, extrap_method)
%TIMETOBOUNDARY Calculates Time-to-Boundary (tau)
%
% ARGUMENTS
% r_x - ML position of the center of pressure or center of mass
%
% r_y - AP position of the center of pressure or center of mass
%
% dt - Unit change in time between samples (1/fs)
%
% bounds - Matrix of boundary coordinates (x in the first column, y in the
% second). Boundary coordinates should be entered in ordered clockwise, but
% need not begin with any particular coordinate. If there are n boundary points (and
% therefore n boundaries) this matrix should be size n x 2.
%
% extrap_method - Determines the method for estimating tau. 
%          The default method is Slobounov (method 2)
%
%     Method 1 (Riccio)
%     Linear Equation:
%     C * tau + D = 0
%     where tau = virtual time-to-boundary
%       C = [v_y(t) - m_bound * v_x(t)]
%       D = [(p_y(t) - y_bound) - m_bound * (p_x(t) - x_bound)]
%     
%     Method 2 (Slobounov)
%     Quadratic Equation:
%     B * tau^2 + C * tau + D = 0
%     where tau = virtual time-to-boundary
%       B = [a_y(t) - m_bound * a_x(t)] / 2
%       C = [v_y(t) - m_bound * v_x(t)]
%       D = [(p_y(t) - y_bound) - m_bound * (p_x(t) - x_bound)]
%     
%     Method 3 (Jerk)
%     Cubic Equation:
%     A * tau^3 + B * tau^2 + C * tau + D = 0
%     where tau = virtual time-to-boundary
%       A = [j_y(t) - m_bound * j_x(t)] / 3
%       B = [a_y(t) - m_bound * a_x(t)] / 2
%       C = [v_y(t) - m_bound * v_x(t)]
%       D = [(p_y(t) - y_bound) - m_bound * (p_x(t) - x_bound)]
%
% RETURNS
% ttb - Time series of minimum Time-to-Boundary in seconds. Minimum TtB is
% determined by taking smallest of the positive, real roots to the
% polynomial equations corresponding to the time it would take the
% extrapolated CoP/CoM trajectory to reach a boundary based on the
% instantaneous position, velocity, and higher-order derivatives.
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

% Argument validation ----------------------------------------------------%

% The first four arguments are required.
if nargin < 4 
    error('Please provide the CoP/CoM coordinates, dt, and boundary positions.');
end

% The default method is 2 = Slobounov
if nargin == 4
    extrap_method = 2;
end

% Compute velocity and higher order derivatives --------------------------%

% Velocity
v_x = gradient(r_x, dt);
v_y = gradient(r_y, dt);

% Acceleration
a_x = gradient(v_x, dt);
a_y = gradient(v_y, dt);

% Jerk
j_x = gradient(a_x, dt);
j_y = gradient(a_y, dt);

% Define slope (s), coordinates (x_bound, y_bound), and polynomial 
% coefficients (A, B, C, D) for each boundary. ---------------------------%

% Length of the data
n_samples = length(r_x);

% Number of boundaries
[n_boundaries, ~] = size(bounds);

% Boundary slopes
slope = nan(n_boundaries, 1);

% Boundary coordinates
x_bound = nan(size(slope));
y_bound = nan(size(slope));

% Polynomial coefficients
coef_A = cell(n_boundaries, 1);
coef_B = cell(n_boundaries, 1);
coef_C = cell(n_boundaries, 1);
coef_D = cell(n_boundaries, 1);

% Polynomial matrix
poly_coef = cell(n_boundaries, 1);

% Roots matrix
poly_roots = cell(n_boundaries, 1);

% For each boundary
for i = 1:n_boundaries
    
    % Compute the boundary slope. For the first boundary, use the first
    % and last boundary positions. For the other boundaries, use the current
    % and previous boundary positions.
    if i == 1
        slope(i) = (bounds(n_boundaries, 2) - bounds(i, 2)) / (bounds(n_boundaries, 1) - bounds(i, 1));
    else
        slope(i) = (bounds(i-1, 2) - bounds(i, 2)) / (bounds(i-1, 1) - bounds(i, 1));
    end

    % If the boundary is vertical (i.e., slope is inf), set the slope to a
    % large number. This should never be needed in practice because it is
    % very unlikely that a boundary would be perfectly vertical.
    if slope(i) == inf
        slope(i) = 1e16;
    end

    % Assign boundary coordinates (x_bound, y_bound). Either of the points 
    % used to define the boundary can be used.
    x_bound(i) = bounds(i, 1);
    y_bound(i) = bounds(i, 2);
    
    % Determine polynomial coefficients (A, B, C, D)
    coef_A{i} = (j_y - slope(i) .* j_x) ./ 3;
    coef_B{i} = (a_y - slope(i) .* a_x) ./ 2;
    coef_C{i} = (v_y - slope(i) .* v_x);
    coef_D{i} = (r_y - y_bound(i)) - slope(i) .* (r_x - x_bound(i));
    
    % Set unnecessary coefficients to 0 based on extrapolation method
    if extrap_method == 1
        coef_A{i} = zeros(n_samples, 1);
        coef_B{i} = zeros(n_samples, 1);
    elseif extrap_method == 2
        coef_A{i} = zeros(n_samples, 1);
    end
    
    % Create polynomial and root matrices
    poly_coef{i} = [coef_A{i} coef_B{i} coef_C{i} coef_D{i}];
    poly_roots{i} = nan(n_samples, extrap_method);
    
end

% Create minimum TtB vector, boundary TtB matrix, boundary crossing vector,
% and matrix for all boundary roots.
ttb = zeros(n_samples, 1);
ttb_bound = zeros(n_samples, n_boundaries);
bound_crossed = zeros(n_samples, 1);
all_roots = zeros(n_samples, n_boundaries * extrap_method);

% Compute TtB ------------------------------------------------------------%

% For each time point
for t = 1:n_samples
      
    % For each boundary
    for i = 1:n_boundaries
        
        % Find roots of polynomial and add to full matrix
        poly_roots{i}(t, :) = roots(poly_coef{i}(t, :))';
        all_roots(t, (extrap_method*i-(extrap_method-1)):extrap_method*i) = poly_roots{i}(t, :);
        
        % Find real roots greater than or equal to 0
        idx_real_pos = find(poly_roots{i}(t, :) >= 0 & imag(poly_roots{i}(t, :)) == 0);
        
        % Find minimum time-to-boundary for this boundary
        if isempty(idx_real_pos)
            ttb_bound(t, i) = NaN;
        else
            ttb_bound(t, i) = min(poly_roots{i}(t, idx_real_pos));
        end
               
        % If TtB to all boundaries has been completed
        if i == n_boundaries
             
            % Find real roots greater than or equal to 0
            idx_real_pos = find(all_roots(t, :) >= 0 & imag(all_roots(t, :)) == 0);
            
            % Find minimum time-to-boundary from the positive real roots
            [ttb(t), min_idx] = min(all_roots(t, idx_real_pos));
            
            % Identify the boundary that was crossed. Divides by the
            % number of roots per boundary and increments to the next 
            % integer to get the corresponding boundary number.
            bound_crossed(t) = ceil(idx_real_pos(min_idx) / extrap_method);
            
        end
    end    
end

% Compute the percent of minimum time-to-boundary to each boundary
bound_percent = ttbBoundaryPercent(bound_crossed, n_boundaries);

end
