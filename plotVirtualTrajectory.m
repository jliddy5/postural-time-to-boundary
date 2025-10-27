function [ ] = plotVirtualTrajectory( p_x, p_y, dt, x_bos, y_bos, ttc, extrapMethod)
% Plots the base of support, center of pressure, and  virtual trajectory.

% Compute higher order derivatives ----------------------------------------%
% Velocity
v_x = cent_diff3(p_x, dt);
v_y = cent_diff3(p_y, dt);

% Acceleration
if extrapMethod == 2 || extrapMethod == 3
    a_x = cent_diff3(v_x, dt);
    a_y = cent_diff3(v_y, dt);
else 
    a_x = zeros(size(v_x));
    a_y = zeros(size(v_y));
end

% Jerk
if extrapMethod == 3
    j_x = cent_diff3(a_x, dt);
    j_y = cent_diff3(a_y, dt);
else 
    j_x = zeros(size(v_x));
    j_y = zeros(size(v_y));
end

% Indices to Plot Extrapolated Trajectory
idx = 950:955;

figure('color','white');
plot(x_bos, y_bos,'k');
xlabel('ML CoP (mm)');
ylabel('AP CoP (mm)');
axis([0 525 0 525]);
hold on;
scatter(x_bos,y_bos,'k','filled');
plot(p_x,p_y,'k');

% For each point in time.
for i = idx

    % Tau time series and virtual trajectory
    tau = 0:0.001:ttc(i);
    vt_x = zeros(size((tau)));
    vt_y = zeros(size(tau));

    % For each tau value.
    for j = 1:length(tau)
        vt_x(j) = p_x(i) + v_x(i) * tau(j) + 0.5 * a_x(i) * (tau(j) ^ 2) + (1/3) * j_x(i) * tau(j) * tau(j) * tau(j);
        vt_y(j) = p_y(i) + v_y(i) * tau(j) + 0.5 * a_y(i) * (tau(j) ^ 2) + (1/3) * j_y(i) * tau(j) * tau(j) * tau(j);
    end

    plot(vt_x, vt_y, 'k');
    hold on;

end

end
