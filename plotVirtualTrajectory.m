function [] = plotVirtualTrajectory(p_x, p_y, dt, x_bos, y_bos, ttb, extrap_method)
%PLOTVIRTUALTRAJECTORY Plots the base of support, center of pressure, and virtual trajectory.
%
% ARGUMENTS
% p_x - ML position of CoP/CoM
% p_y - AP position of CoP/CoM
% dt - Time step
% x_bos - Base of support x-coordinates (closed polygon)
% y_bos - Base of support y-coordinates (closed polygon)
% ttb - Time-to-boundary values
% extrap_method - Extrapolation method (1=Riccio, 2=Slobounov, 3=Jerk)

% Compute higher order derivatives ----------------------------------------%

% Velocity
v_x = gradient(p_x, dt);
v_y = gradient(p_y, dt);

% Acceleration
if extrap_method == 2 || extrap_method == 3
    a_x = gradient(v_x, dt);
    a_y = gradient(v_y, dt);
else 
    a_x = zeros(size(v_x));
    a_y = zeros(size(v_y));
end

% Jerk
if extrap_method == 3
    j_x = gradient(a_x, dt);
    j_y = gradient(a_y, dt);
else 
    j_x = zeros(size(v_x));
    j_y = zeros(size(v_y));
end

% Indices to plot extrapolated trajectory
plot_indices = 950:955;

f = figure();
theme(f, 'light');
plot(x_bos, y_bos,'k');
xlabel('ML CoP (mm)');
ylabel('AP CoP (mm)');
axis([0 525 0 525]);
hold on;
scatter(x_bos, y_bos, 'k', 'filled');
plot(p_x, p_y, 'k');

% For each point in time
for i = plot_indices

    % Tau time series and virtual trajectory
    tau = 0:0.001:ttb(i);
    vt_x = zeros(size((tau)));
    vt_y = zeros(size(tau));

    % For each tau value
    for j = 1:length(tau)
        vt_x(j) = p_x(i) + v_x(i) * tau(j) + 0.5 * a_x(i) * (tau(j)^2) + (1/3) * j_x(i) * tau(j)^3;
        vt_y(j) = p_y(i) + v_y(i) * tau(j) + 0.5 * a_y(i) * (tau(j)^2) + (1/3) * j_y(i) * tau(j)^3;
    end

    plot(vt_x, vt_y, 'k');
    hold on;

end

end
