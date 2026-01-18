clc; clear all; close all

%% Read in saved states
% load('states_1day_10el.mat')
% load('states_1day_20d2el.mat')
% load('states_3day_20d2el.mat')
% load('states_1day_20d2el_wd.mat')
load('states_3day_20d2el_wd.mat')



%% Initialize Minimum LOS distance
num_sats = 3;
num_gs = 3;

% initialize satellites
sats_lunar_pos = zeros(3,length(times),num_sats);
sats_lunar_pos(:,:,1) = sat_lunar1_pos;
sats_lunar_pos(:,:,2) = sat_lunar2_pos;
sats_lunar_pos(:,:,3) = sat_lunar3_pos;

% initialize ground stations
gss_pos = zeros(3,length(times),num_gs);
gss_pos(:,:,1) = gs1_pos;
gss_pos(:,:,2) = gs2_pos;
gss_pos(:,:,3) = gs3_pos;

% initialize all access objects
acs = [ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9];

% find minimum LOS distance
min_LOS_dist_arr = zeros([1,length(times_datetime)]);
acs_index = zeros([1,length(times)]);
num_LOS_pairs = zeros([1,length(times_datetime)]);
for i=1:length(times_datetime)

    % check if LOS to gs
    isLOS = accessStatus(acs,times_datetime(i)); 

    min_LOS_dist = Inf;
    counter = 1;
    for j=1:num_gs
        for k=1:num_sats
            if isLOS(counter) && is_lunar_sat_accessible(sats_lunar_pos(:,i,k),gss_pos(:,i,j),moon_pos(:,i),90)
                distance = norm(sats_lunar_pos(:,i,k) - gss_pos(:,i,j));
                num_LOS_pairs(i) = num_LOS_pairs(i) + 1;
                if distance < min_LOS_dist
                    min_LOS_dist = distance;
                    acs_index(i) = counter;
                end
            end
            counter = counter + 1;
        end
    end
   min_LOS_dist_arr(i) =  min_LOS_dist;
end


%% Initialize Lunar Satellite Cone Angles and GS Elevation Angles

lunar_sat_cone_angles = zeros([1,length(times_datetime)]);
gs_el_angles = zeros([1,length(times_datetime)]);
for i=1:length(times)

    if acs_index(i) == 0
        lunar_sat_cone_angles(i) = NaN;
        gs_el_angles(i) = NaN;
    else

        if mod(acs_index(i),3) == 0
            sc_interest = 3;
        else
            sc_interest = mod(acs_index(i),3);
        end
        gs_interest = ceil(acs_index(i)/3);

        % initialize lunar cone angle
        lunar_sat_zenith_vector = sats_lunar_pos(:,i,sc_interest) - moon_pos(:,i);
        lunar_sat_to_gs_vector = gss_pos(:,i,gs_interest) - sats_lunar_pos(:,i,sc_interest);
        lunar_sat_cone_angles(i) = angle_between_vectors(lunar_sat_zenith_vector,lunar_sat_to_gs_vector);

        % initialize gs elevation angle
        gs_zenith_vector = gss_pos(:,i,gs_interest);    % center of Earth is at (0,0,0)
        gs_to_lunar_sat_vector = sats_lunar_pos(:,i,sc_interest) - gss_pos(:,i,gs_interest);
        gs_el_angles(i) = 90 - angle_between_vectors(gs_zenith_vector,gs_to_lunar_sat_vector);
    end

end


%% Compute Link Budget

% constants
c = 3e8;
k = 1.38e-23;

rss = zeros([1,length(times)]);
snr = zeros([1,length(times)]);
R = zeros([1,length(times)]);
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
    pointing_weather_loss = 4;                    % worst case loss with pointing and weather

    % receiver segment
    rss(i) = Pt - tx_cable_loss + Gt - FSPL - pointing_weather_loss;    % dBm
    Gr = 79-0.3;                            % worst case from DSN user guide
    rx_cable_loss = 0;                      % dB (DSN waveguides have high performance)
    Ps = rss(i) + Gr  - rx_cable_loss;      % dBm

    % noise analysis
    T_sys = 44;                             % K, from DSN user guide, 40K + 10% worst case
    Pnoise_linear = k*T_sys*bw;             % W
    Pnoise = 10*log10(Pnoise_linear) + 30;  % dBm
    snr(i) = Ps - Pnoise;

    % achievable bitrate
    snr_linear = 10^(snr(i)/10);
    R(i) = bw*log2(1+snr_linear);

end


%% Plots

% minimum LOS distance plot
figure
plot(times/3600,(min_LOS_dist_arr/1e3)/R_earth,LineWidth=2)
xlim([0,max(times/3600)])
xlabel('Time [hrs]')
ylabel('Distance [Earth Radii]')
% title('Minimum LOS Distance vs Time')
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
plot(times/3600,lunar_sat_cone_angles,LineWidth=2,DisplayName='Lunar Satellite Cone Angle')
hold on
plot(times/3600,gs_el_angles,LineWidth=2,DisplayName='Ground Station El Angle')
xlim([0,max(times/3600)])
xlabel('Time [hrs]')
ylabel('Angle [deg]')
% title('Cone Angle and Elevation Angle vs Time')
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

% SNR plot
figure
hold on
plot(times/3600,snr,LineWidth=2)
xlim([0,max(times/3600)])
xlabel('Time [hrs]')
ylabel('SNR [dB]')
% title('SNR vs Time')
grid on
hold off

% datarate plot
M_OPQSK = 4;
M_GMSK = 2;
coding_rates = [1/2,2/3,4/5,7/8];
rate_limits_OPQSK = zeros([1,length(coding_rates)]);
rate_limits_GMSK = zeros([1,length(coding_rates)]);

for i=1:length(coding_rates)
    rate_limits_OPQSK(i) = bw*coding_rates(i)*log2(M_OPQSK);
    rate_limits_GMSK(i) = bw*coding_rates(i)*log2(M_GMSK);
end

figure
semilogy(times/3600,R/1e6,LineWidth=2)
xlim([0,max(times/3600)])
xlabel('Time [hrs]')
ylabel('Datarate [Mbps]')
% title('Maximum Datarate (Shannon Limit)')
grid on
hold off