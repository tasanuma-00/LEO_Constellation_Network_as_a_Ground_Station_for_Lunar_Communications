clc; clear all; close all

%% Read in saved states
load('states_1day_10el.mat')


%% Plot in Earth Reference Frame

% Initialize the figure and 3D plot
figure;
view(3);
hold on;
grid on;
xlim([-2*d_earth_moon, 2*d_earth_moon]);
ylim([-2*d_earth_moon, 2*d_earth_moon]);
zlim([-2*d_earth_moon, 2*d_earth_moon]);

% Pre-create the plot object for efficiency
h = plot3(NaN, NaN, NaN, 'b.', 'MarkerSize', 10);
h2 = plot3(NaN, NaN, NaN, 'r.', 'MarkerSize', 10);
h3 = plot3(NaN, NaN, NaN, 'g.', 'MarkerSize', 10);

% Initialize data storage for the plot
xyz_data_sat_lunar1 = sat_lunar1_pos(:,1);
xyz_data_moon = moon_pos(:,1);
xyz_data_gs1 = gs1_pos(:,1);

% Loop to update the plot
for i = 2:length(times)

    % Append the new point to the data arrays
    xyz_data_sat_lunar1 = [[xyz_data_sat_lunar1(1,:), sat_lunar1_pos(1, i) / 1e3]; ...
        [xyz_data_sat_lunar1(2,:), sat_lunar1_pos(2, i) / 1e3]; ...
        [xyz_data_sat_lunar1(3,:), sat_lunar1_pos(3, i) / 1e3]];

    xyz_data_moon = [[xyz_data_moon(1,:), moon_pos(1, i) / 1e3]; ...
        [xyz_data_moon(2,:), moon_pos(2, i) / 1e3]; ...
        [xyz_data_moon(3,:), moon_pos(3, i) / 1e3]];
    
    xyz_data_gs1 = [[xyz_data_gs1(1,:), gs1_pos(1, i) / 1e3]; ...
        [xyz_data_gs1(2,:), gs1_pos(2, i) / 1e3]; ...
        [xyz_data_gs1(3,:), gs1_pos(3, i) / 1e3]];

    
    % Update the data of the plot object
    num_trail = 10000;
    if i > num_trail
        set(h, 'XData', xyz_data_sat_lunar1(1,end-num_trail:end), ...
            'YData', xyz_data_sat_lunar1(2,end-num_trail:end), ...
            'ZData', xyz_data_sat_lunar1(3,end-num_trail:end));
        set(h2, 'XData', xyz_data_moon(1,end-num_trail:end), ...
            'YData', xyz_data_moon(2,end-num_trail:end), ...
            'ZData', xyz_data_moon(3,end-num_trail:end));
        set(h3, 'XData', xyz_data_gs1(1,end-num_trail:end), ...
            'YData', xyz_data_gs1(2,end-num_trail:end), ...
            'ZData', xyz_data_gs1(3,end-num_trail:end));

    else
        set(h, 'XData', xyz_data_sat_lunar1(1,:), ...
            'YData', xyz_data_sat_lunar1(2,:), ...
            'ZData', xyz_data_sat_lunar1(3,:));
        set(h2, 'XData', xyz_data_moon(1,:), ...
            'YData', xyz_data_moon(2,:), ...
            'ZData', xyz_data_moon(3,:));
        set(h3, 'XData', xyz_data_gs1(1,:), ...
            'YData', xyz_data_gs1(2,:), ...
            'ZData', xyz_data_gs1(3,:));
    end

    % Refresh the plot
    drawnow limitrate; % Efficiently updates the figure
end
hold off;


%% Plot in Moon Reference Frame

% Initialize the figure and 3D plot
figure;
view(3);
hold on;
grid on;

% Set axis limits
xlim([-2*R_earth, 2*R_earth]);
ylim([-2*R_earth, 2*R_earth]);
zlim([-2*R_earth, 2*R_earth]);

% Pre-create the plot object for efficiency
h = plot3(NaN, NaN, NaN, 'b.', 'MarkerSize', 10);
h2 = plot3(NaN, NaN, NaN, 'r.', 'MarkerSize', 10);
h3 = plot3(NaN, NaN, NaN, 'g.', 'MarkerSize', 10);


% Initialize data storage for the plot
x_data_sat_lunar1 = [];
y_data_sat_lunar1 = [];
z_data_sat_lunar1 = [];

x_data_sat_lunar2 = [];
y_data_sat_lunar2 = [];
z_data_sat_lunar2 = [];

x_data_sat_lunar3 = [];
y_data_sat_lunar3 = [];
z_data_sat_lunar3 = [];

x_data_moon = [];
y_data_moon = [];
z_data_moon = [];


% Loop to update the plot
for i = 1:length(times)
    % Append the new point to the data arrays
    x_data_moon = [x_data_moon, moon_pos(1, i) / 1e3];
    y_data_moon = [y_data_moon, moon_pos(2, i) / 1e3];
    z_data_moon = [z_data_moon, moon_pos(3, i) / 1e3];

    x_data_sat_lunar1 = [x_data_sat_lunar1, sat_lunar1_pos(1, i) / 1e3 - x_data_moon(i)];
    y_data_sat_lunar1 = [y_data_sat_lunar1, sat_lunar1_pos(2, i) / 1e3 - y_data_moon(i)];
    z_data_sat_lunar1 = [z_data_sat_lunar1, sat_lunar1_pos(3, i) / 1e3 - z_data_moon(i)];

    x_data_sat_lunar2 = [x_data_sat_lunar2, sat_lunar2_pos(1, i) / 1e3 - x_data_moon(i)];
    y_data_sat_lunar2 = [y_data_sat_lunar2, sat_lunar2_pos(2, i) / 1e3 - y_data_moon(i)];
    z_data_sat_lunar2 = [z_data_sat_lunar2, sat_lunar2_pos(3, i) / 1e3 - z_data_moon(i)];

    x_data_sat_lunar3 = [x_data_sat_lunar3, sat_lunar3_pos(1, i) / 1e3 - x_data_moon(i)];
    y_data_sat_lunar3 = [y_data_sat_lunar3, sat_lunar3_pos(2, i) / 1e3 - y_data_moon(i)];
    z_data_sat_lunar3 = [z_data_sat_lunar3, sat_lunar3_pos(3, i) / 1e3 - z_data_moon(i)];
    
    num_trail = 10;
    if i > 100
        set(h, 'XData', x_data_sat_lunar1(end-num_trail:end), ...
                'YData', y_data_sat_lunar1(end-num_trail:end), ...
                'ZData', z_data_sat_lunar1(end-num_trail:end));
        set(h2, 'XData', x_data_sat_lunar2(end-num_trail:end), ...
                'YData', y_data_sat_lunar2(end-num_trail:end), ...
                'ZData', z_data_sat_lunar2(end-num_trail:end));
        set(h3, 'XData', x_data_sat_lunar3(end-num_trail:end), ...
                'YData', y_data_sat_lunar3(end-num_trail:end), ...
                'ZData', z_data_sat_lunar3(end-num_trail:end));

    else
        set(h, 'XData', x_data_sat_lunar1, ...
                'YData', y_data_sat_lunar1, ...
                'ZData', z_data_sat_lunar1);
        set(h2, 'XData', x_data_sat_lunar2, ...
                'YData', y_data_sat_lunar2, ...
                'ZData', z_data_sat_lunar2);
        set(h3, 'XData', x_data_sat_lunar3, ...
                'YData', y_data_sat_lunar3, ...
                'ZData', z_data_sat_lunar3);
    end
    
    % Refresh the plot
    drawnow limitrate; % Efficiently updates the figure
    pause(1e-3)
end
hold off;
