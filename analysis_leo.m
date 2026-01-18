% Create a Satellite Scenario
clc; clear; close all;

%%
% load('states_1day_20d2el_wd.mat')
load('states_3day_20d2el_wd.mat')

%% Initialize Minimum LOS distance
num_lunar_sats = 3;
num_leo_sats = 24;
earth_pos = [0;0;0];

% initialize satellites
sats_lunar_pos = zeros(3,length(times),num_lunar_sats);
sats_lunar_pos(:,:,1) = sat_lunar1_pos;
sats_lunar_pos(:,:,2) = sat_lunar2_pos;
sats_lunar_pos(:,:,3) = sat_lunar3_pos;

% find minimum LOS distance
min_LOS_dist_arr = zeros([1,length(times_datetime)]);
indices_arr = zeros(2,length(times));
num_LOS_pairs = zeros([1,length(times_datetime)]);
for i=1:length(times_datetime)
    

    min_LOS_dist = Inf;
    counter = 1;
    for j=1:num_leo_sats
        for k=1:num_lunar_sats

            % check if LOS
            is_lunar_LOS = is_lunar_sat_accessible(sats_lunar_pos(:,i,k),leos_pos(:,i,j),moon_pos(:,i),90);
            is_leo_LOS = is_leo_sat_accessible(leos_pos(:,i,j),sats_lunar_pos(:,i,k),earth_pos,90);

            % initialize minimum distance LOS
            if is_lunar_LOS && is_leo_LOS
                distance = norm(sats_lunar_pos(:,i,k) - leos_pos(:,i,j));
                num_LOS_pairs(i) = num_LOS_pairs(i) + 1;
                if distance < min_LOS_dist
                    min_LOS_dist = distance;
                    indices_arr(:,i) = [j,k];
                end
            end
        end
    end
   min_LOS_dist_arr(i) =  min_LOS_dist;
end




%% Initialize Lunar Satellite Cone Angles and GS Elevation Angles

lunar_sat_cone_angles = zeros([1,length(times_datetime)]);
leo_sat_cone_angles = zeros([1,length(times_datetime)]);
for i=1:length(times)

    if indices_arr(1,i) == 0
        lunar_sat_cone_angles(i) = NaN;
        leo_sat_cone_angles(i) = NaN;
    else

        % indices of lunar sat and leo sat pair
        j = indices_arr(1,i);
        k = indices_arr(2,i);

        % initialize lunar cone angle
        lunar_sat_zenith_vector = sats_lunar_pos(:,i,k) - moon_pos(:,i);
        lunar_sat_to_leo_sat_vector = leos_pos(:,i,j) - sats_lunar_pos(:,i,k);
        lunar_sat_cone_angles(i) = angle_between_vectors(lunar_sat_zenith_vector,lunar_sat_to_leo_sat_vector);

        % initialize gs elevation angle
        leo_sat_zenith_vector = leos_pos(:,i,j);    % center of Earth is at (0,0,0)
        leo_sat_to_lunar_sat_vector = sats_lunar_pos(:,i,k) - leos_pos(:,i,j);
        leo_sat_cone_angles(i) = angle_between_vectors(leo_sat_zenith_vector,leo_sat_to_lunar_sat_vector);

    end

end


%% Compute Link Budget

% constants
c = 3e8;
k = 1.38e-23;

rss = zeros([1,length(times)]);
for i=1:length(times)

    % comm parameters
    freq = 22.85e9;
    bw = 0.6e9;

    % transmitter segment
    Pt_linear = 0.3;                % W
    Pt = 10*log10(Pt_linear) + 30;  % dBm
    tx_cable_loss = 0.25;           % dB
    d_antenna = 1;                  % m
    antenna_efficiency = 0.8;       
    Gt_linear = antenna_efficiency * (pi*d_antenna*freq/c)^2;   % assume parabolic antenna
    Gt = 10*log10(Gt_linear);       % dB

    % propagation segment
    FSPL = fspl(min_LOS_dist_arr(i), c/freq);     % dB

    % receiver segment
    rss(i) = Pt - tx_cable_loss + Gt - FSPL;        % dBm

end

%% Plots

% minimum LOS distance plot
figure
plot(times/3600,(min_LOS_dist_arr/1e3)/R_earth,LineWidth=2)
xlim([0,max(times/3600)])
xlabel('Time [hrs]')
ylabel('Distance [Earth Radii]')
% title('Minimum LOS Distance')
grid on

% num LOS pairs plot
figure
plot(times/3600,num_LOS_pairs,LineWidth=2)
xlim([0,max(times/3600)])
xlabel('Time [hrs]')
ylabel('Number of LOS Links')
% title('Number of Possible LOS Links vs Time')
grid on

% cone angle plot
figure
plot(times/3600,lunar_sat_cone_angles,LineWidth=2,DisplayName='Lunar Satellite')
hold on
plot(times/3600,leo_sat_cone_angles,LineWidth=2,DisplayName='LEO Satellite')
xlim([0,max(times/3600)])
xlabel('Time [hrs]')
ylabel('Cone Angle [deg]')
% title('Cone Angles vs Time')
legend(Location='best')
grid on

% RSS plot
figure
hold on
plot(times/3600,rss,LineWidth=2)
xlim([0,max(times/3600)])
xlabel('Time [hrs]')
ylabel('Incident Signal Power [dB]')
% title('Incident Signal Power vs Time')
grid on
hold off