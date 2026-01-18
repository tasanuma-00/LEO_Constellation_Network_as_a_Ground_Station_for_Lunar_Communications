clc; clear all; close all

%%
function r_eci = get_gs_states(lat,lon,alt,utc)
% converts from geodetic -> ecef -> eci
% returns coordinates in eci (meters)

% Convert geodetic coordinates (lat, lon, alt) to ECEF (X, Y, Z)
wgs84 = wgs84Ellipsoid;
[x,y,z] = geodetic2ecef(wgs84,lat,lon,alt);
r_ecef = [x,y,z];
v_ecef = [0,0,0];

% Convert ECEF coordinates to ECI
r_eci = ecef2eci(utc, r_ecef, v_ecef);

end

%% Constants

% Define geodetic coordinates (latitude, longitude, altitude) in degrees and meters
gs_lat = 35.426667; % Latitude of ground station in degrees
gs_lon = -116.89;   % Longitude of ground station in degrees
gs_alt = 0;         % Altitude above sea level in km

% UTC time (for Earth rotation calculation)
utc_time = datetime(2024, 11, 11, 0, 0, 0, 'TimeZone', 'UTC');

dt = 60;                    % Time step in seconds
total_time = 3*24*3600;     % Total simulation time in seconds (3 days)
times = 0:dt:total_time;
start_time = datetime('now');
stop_time = start_time + seconds(total_time);
times_datetime = start_time:seconds(60):stop_time;

r_eci = zeros([3,length(times_datetime)]);
for i=1:length(times_datetime)
    r_eci(:,i) = get_gs_states(gs_lat,gs_lon,gs_alt,times_datetime(i));
end

