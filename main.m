clc; clear all; close all

%% Constants
mu_moon = 4902.8;       % Moon's gravitational parameter in km^3/s^2
R_moon = 1737.4;        % Moon radius in km
R_earth = 6371;         % Earth radius in km
d_earth_moon = 384400;  % Earth-Moon distance in km


%% Construct Satellite Scenario

% initialize time
dt = 60;                    % Time step in seconds
total_time = 3*24*3600;     % Total simulation time in seconds (3 days)
times = 0:dt:total_time;
start_time = datetime('2024-01-01','InputFormat','yyyy-MM-dd');
stop_time = start_time + seconds(total_time);
times_datetime = start_time:seconds(60):stop_time;

% Initialize satelliteScenario with numerical propagator
sc = satelliteScenario(start_time,stop_time,dt);
numericalPropagator(sc, ...
    ODESet=odeset(RelTol=1e-6,AbsTol=1e-6,MaxStep=60), ...
    IncludeThirdBodyGravity=true, ...
    ThirdBodyGravitySource=["Sun" "Mercury" "Venus" "Moon" "Mars" "Jupiter" "Saturn" "Uranus" "Neptune" "Pluto"], ...
    GravitationalPotentialModel="point-mass")


%% Moon

% construct time table for moon using JPL ephemeris
leapSeconds = 37;    % s
ttMinusTAI = 32.184; % s
terrestrialTime = times_datetime + seconds(leapSeconds + ttMinusTAI);
tdbJD = tdbjuliandate([ ...
    terrestrialTime.Year' ...
    terrestrialTime.Month' ...
    terrestrialTime.Day' ...
    terrestrialTime.Hour' ...
    terrestrialTime.Minute' ...
    terrestrialTime.Second']);
pMoonkm = planetEphemeris(tdbJD,"earth","moon"); % km
pMoon = convlength(pMoonkm,'km','m');            % m
pMoonTT = timetable(times_datetime',pMoon);

% initialize moon as a satellite object
moon = satellite(sc,pMoonTT,Name="Moon");
moon.Orbit.LineColor="red";
moon.MarkerColor="red";


%% Create Lunar Satellites

% Orbit Parameters w.r.t Moon
a = 6142.4;             % Semi-major axis in km
e = 0;                  % Eccentricity
i = 0;                  % Inclination in radians
raan = deg2rad(0);      % Right ascension of the ascending node in radians
aop = deg2rad(315);     % Argument of periapsis (not needed for circular orbit)
true_anom = deg2rad(0); % True anomaly in radians
[a, ecc, incl, RAAN, argp, nu] = init_lunar_satellite(a,e,i,raan,aop,true_anom,tdbJD);
sat_lunar1 = satellite(sc, a, ecc, incl, RAAN, argp, nu, ...
    OrbitPropagator="numerical", Name="s/c 1");


% Orbit Parameters w.r.t Moon
a = 6142.4;             % Semi-major axis in km
e = 0.59999;            % Eccentricity
i = deg2rad(57.7);      % Inclination in radians
raan = deg2rad(270);    % Right ascension of the ascending node in radians
aop = deg2rad(270);     % Argument of periapsis (not needed for circular orbit)
true_anom = deg2rad(0); % True anomaly in radians

[a, ecc, incl, RAAN, argp, nu] = init_lunar_satellite(a,e,i,raan,aop,true_anom,tdbJD);
sat_lunar2 = satellite(sc, a, ecc, incl, RAAN, argp, nu, ...
    OrbitPropagator="numerical", Name="s/c 2");

% Orbit Parameters w.r.t Moon
a = 6142.4;             % Semi-major axis in km
e = 0.59999;            % Eccentricity
i = deg2rad(57.7);      % Inclination in radians
raan = deg2rad(0);      % Right ascension of the ascending node in radians
aop = deg2rad(90);      % Argument of periapsis (not needed for circular orbit)
true_anom = deg2rad(0); % True anomaly in radians

[a, ecc, incl, RAAN, argp, nu] = init_lunar_satellite(a,e,i,raan,aop,true_anom,tdbJD);
sat_lunar3 = satellite(sc, a, ecc, incl, RAAN, argp, nu, ...
    OrbitPropagator="numerical", Name="s/c 3");


%% Define ground station on Earth using geodetic coordinates (latitude, longitude, altitude)

% Goldstone
gs_lat1 = 35.426667; % Latitude of ground station in degrees
gs_lon1 = -116.89;   % Longitude of ground station in degrees
gs_alt1 = 0;         % Altitude above sea level in km

% Madrid
gs_lat2 = 40.431261988480294; 
gs_lon2 = -4.247946445841926;
gs_alt2 = 0;

% Canberra
gs_lat3 = -35.402343860700995; 
gs_lon3 = 148.98290398120525;
gs_alt3 = 0;

gs1 = groundStation(sc, gs_lat1, gs_lon1, 'Altitude', gs_alt1*1e3,'MinElevationAngle',20.2);
gs2 = groundStation(sc, gs_lat2, gs_lon2, 'Altitude', gs_alt2*1e3,'MinElevationAngle',20.2);
gs3 = groundStation(sc, gs_lat3, gs_lon3, 'Altitude', gs_alt3*1e3,'MinElevationAngle',20.2);


%% Create LEO Constellation

% Define Walker Delta constellation parameters
numSatellites = 24;    % Total number of satellites
numPlanes = 6;         % Number of orbital planes
phasingFactor = 1;     % Phasing factor
altitude = 600e3;      % Altitude (meters)
inclination = 97.4;    % Inclination (degrees, Sun-synchronous)

% Generate Walker Delta constellation using the deltawalker function
walkerConstellation = walkerDelta(sc, R_earth*1e3+altitude, inclination, ...
                                   numSatellites, numPlanes, phasingFactor);


%% Visualization (doesn't work for complex sims)

% v = satelliteScenarioViewer(sc, ...
%     CameraReferenceFrame="Inertial", ...
%     PlaybackSpeedMultiplier=6000);


%% Obtain States

% moon 
[moon_pos,moon_vel] = states(moon);

% lunar satellites
[sat_lunar1_pos,sat_lunar1_vel] = states(sat_lunar1);
[sat_lunar2_pos,sat_lunar2_vel] = states(sat_lunar2);
[sat_lunar3_pos,sat_lunar3_vel] = states(sat_lunar3);

% LEO satellites
leos_pos = states(walkerConstellation);

% ground stations
gs1_pos = zeros([3,length(times_datetime)]);
gs2_pos = zeros([3,length(times_datetime)]);
gs3_pos = zeros([3,length(times_datetime)]);
for i=1:length(times_datetime)
    gs1_pos(:,i) = get_gs_states(gs_lat1,gs_lon1,gs_alt1,times_datetime(i));
    gs2_pos(:,i) = get_gs_states(gs_lat2,gs_lon2,gs_alt2,times_datetime(i));
    gs3_pos(:,i) = get_gs_states(gs_lat3,gs_lon3,gs_alt3,times_datetime(i));
end


%% Access
% access only accounting for Earth's blockage
ac1 = access(sat_lunar1,gs1);
ac2 = access(sat_lunar2,gs1);
ac3 = access(sat_lunar3,gs1);
ac4 = access(sat_lunar1,gs2);
ac5 = access(sat_lunar2,gs2);
ac6 = access(sat_lunar3,gs2);
ac7 = access(sat_lunar1,gs3);
ac8 = access(sat_lunar2,gs3);
ac9 = access(sat_lunar3,gs3);

intvls1 = accessIntervals(ac1);
intvls2 = accessIntervals(ac2);
intvls3 = accessIntervals(ac3);
intvls4 = accessIntervals(ac4);
intvls5 = accessIntervals(ac5);
intvls6 = accessIntervals(ac6);
intvls7 = accessIntervals(ac7);
intvls8 = accessIntervals(ac8);
intvls9 = accessIntervals(ac9);