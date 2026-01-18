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