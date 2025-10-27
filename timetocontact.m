function [ ttc, ttc_bound, bound_crossed, bound_percent ] = timetocontact( r_x, r_y, dt, bounds, extrapMethod )
%TTC Calculates Time-to-Contact (tau)
%
% ARGUMENTS
% p_x - ML position of the center of pressure or center of mass
%
% p_y - AP position of the center of pressure or center of mass
%
% dt - Unit change in time between samples (1/fs)
%
% bounds - Matrix of boundary coordinates (x in the first column, y in the
% second). Boundary coordinates should be entered in ordered clockwise, but
% need not begin with any particular coordinate. If there are n boundary points (and
% therefore n boundaries) this matrix should be size n x 2.
%
% method - Determines the method for estimating tau. 
%          The default method is Slobounov (method 2)
%
%     Method 1 (Riccio)
%     Linear Equation:
%     C * tau + D = 0
%     where tau = virtual time-to-contact
%       C = [v_y(t) - m_bound * v_x(t)]
%       D = [(p_y(t) - y_bound) - m_bound * (p_x(t) - x_bound)]
%     
%     Method 2 (Slobounov)
%     Quadratic Equation:
%     B * tau^2 + C * tau + D = 0
%     where tau = virtual time-to-contact
%       B = [a_y(t) - m_bound * a_x(t)] / 2
%       C = [v_y(t) - m_bound * v_x(t)]
%       D = [(p_y(t) - y_bound) - m_bound * (p_x(t) - x_bound)]
%     
%     Method 3 (Jerk)
%     Quadratic Equation:
%     A * tau^3 + B * tau^2 + C * tau + D = 0
%     where tau = virtual time-to-contact
%       A = [j_y(t) - m_bound * j_x(t)] / 3
%       B = [a_y(t) - m_bound * a_x(t)] / 2
%       C = [v_y(t) - m_bound * v_x(t)]
%       D = [(p_y(t) - y_bound) - m_bound * (p_x(t) - x_bound)]
%
%RETURNS
% ttc - Time series of minimum Time-to-Contact in seconds. Minimum TtC is
% determined by taking smallest of the positive, real roots to the
% polynomial equations corresponding to the time it would take the
% extrapolated CoP/CoM trajectory to reach a boundary based on the
% instanteous position, velocity, and higher-order derivatives.
%
% bound_ttc - Minimum Time-to-Contact values for each boundary. The minimum
% of the positive, real roots for each point in time is selected. If no
% positive, real roots exist a value of NaN is returned.
%
% bound_perc - Calculates the percentage of virtual contacts to each
% boundary. Uses the ttc_bound_perc function.
%
% ========================================================================%

% Argument validation ----------------------------------------------------%

% The first four arguments are required.
if nargin < 4 
    disp('Error: Please provide the CoP/CoM coordinates, dt, and boundary positions.');
    return;
end

% The default method is 2 = Slobounov
if nargin == 4
    extrapMethod = 2;
end

% Compute velocity and higher order derivatives --------------------------%

% Velocity
v_x = cent_diff3(r_x, dt);
v_y = cent_diff3(r_y, dt);

% Acceleration
a_x = cent_diff3(v_x, dt);
a_y = cent_diff3(v_y, dt);

% Jerk
j_x = cent_diff3(a_x, dt);
j_y = cent_diff3(a_y, dt);

% Define slope (s), coordinates (x_b, y_b), and polynomial coefficients 
% (A, B, C, D) for each boundary. ----------------------------------------%

% Length of the data
N = length(r_x);

% Number of boundaries
[nBoundaries, ~] = size(bounds);

% Boundary slopes
s = nan(nBoundaries,1);

% Boundary coordinates
x_b = nan(size(s));
y_b = nan(size(s));

% Polynomial coefficients
A = cell(nBoundaries,1);
B = cell(nBoundaries,1);
C = cell(nBoundaries,1);
D = cell(nBoundaries,1);

% Polynomial matrix
p = cell(nBoundaries,1);

% Roots matrix
r = cell(nBoundaries,1);

% For each boundary
for i = 1:nBoundaries
    
    % Compute the boundary slope (s). For the first boundary, use the first
    % and last boundary positions. For the other boundaries, use the current
    % and previous boundary positions.
    if i == 1
        s(i) = (bounds(nBoundaries,2) - bounds(i,2)) / ( bounds(nBoundaries,1) - bounds(i,1));
    else
        s(i) = (bounds(i-1,2) - bounds(i,2)) / ( bounds(i-1,1) - bounds(i,1));
    end

    % If the boundary is vertical (i.e., s is inf), set the slope to a
    % large number. This should never be needed in practice because it is
    % very unlikely that a boundary would be perfectly vertical.
    if s(i) == inf
        s(i) = 1e16;
    end

    % Assign boundary coordinates (x_b, y_b). Either of the points used to
    % define the boundary can be used.
    x_b(i) = bounds(i,1);
    y_b(i) = bounds(i,2);
    
    % Determine polynomial coefficients (A, B, C, D)
    A{i} = (j_y - s(i).* j_x) ./ 3;
    B{i} = (a_y - s(i).* a_x) ./ 2;
    C{i} = (v_y - s(i).* v_x);
    D{i} = (r_y - y_b(i)) - s(i) .* (r_x - x_b(i));
    
    % Set unnecessary coefficients to 0
    if extrapMethod == 1
        A{i} = zeros(N,1);
        B{i} = zeros(N,1);
    elseif extrapMethod == 2
        A{i} = zeros(N,1);
    end
    
    % Create polynomial and root matrices.
    p{i} = [A{i} B{i} C{i} D{i}];
    r{i} = nan(N,extrapMethod);
    
end %end boundary

% Create minimum TtC vector, boundary TtC matrix, boundary crossing vector,
% and matrix for all boundary roots.
ttc = zeros(N,1);
ttc_bound = zeros(N,nBoundaries);
bound_crossed = zeros(N,1);
rts = zeros(N, nBoundaries*extrapMethod);

% Compute TtC ------------------------------------------------------------%

% For each time point
for t = 1:N
      
    % For each boundary
    for i = 1:nBoundaries
        
        % Find roots of polynomial and add to full matrix
        r{i}(t,:) = roots(p{i}(t,:))';
        rts(t, (extrapMethod*i-(extrapMethod-1)):extrapMethod*i ) =  r{i}(t,:);
        
        % Find real roots greater than or equal to 0. 
        idx = find(r{i}(t,:) >= 0 & imag(r{i}(t,:)) == 0);
        
        % Find minimum time to contact for boundary
        if isempty(idx)
            ttc_bound(t,i) = NaN;
        else
            ttc_bound(t,i) = min(r{i}(t,idx));
        end
               
        % If TtC to all boundaries has been completed.
        if i == nBoundaries
             
            % Find real roots greater than or equal to 0. 
            idx = find(rts(t,:) >= 0 & imag(rts(t,:)) == 0 );
            
            % Find minimum time to contact from the positive real roots.
            [ ttc(t), min_idx ] = min(rts(t,idx));
            
            % Identify the boundary that was crossed. Divides by the
            % number of roots per boundary and increments to the next 
            % integer to get the corresponding boundary number.
            bound_crossed(t) = ceil(idx(min_idx)/extrapMethod);
            
        end
    end    
end

% Compute the percent of minimum time-to-contact to each boundary.
bound_percent = ttcBoundaryPercent( bound_crossed, nBoundaries);

end