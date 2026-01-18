function [a, ecc, incl, RAAN, argp, nu] = init_lunar_satellite(a,e,i,raan,aop,true_anom,tdbJD)
    
    % gravitational parameter of moon
    mu_moon = 4902.8; % Moon's gravitational parameter in km^3/s^2

    % Calculated orbit parameters
    p = a*(1-e^2);
    true_long = true_anom + aop + raan;
    arg_lat = true_anom + aop;
    long_per = aop + raan;

    % Calculate position and velocity vectors in the Moon-centered inertial frame
    [r_sat_rel_moon, v_sat_rel_moon] = orb2rv(p,e,i,raan,aop,true_anom,true_long,arg_lat,long_per,mu_moon);  % km,km/s

    % Obtain the position and velocity of the Moon relative to Earth
    [r_moon, v_moon] = planetEphemeris(tdbJD(1), "earth", "moon");  % km, km/s

    % Transform to Earth-centered frame
    r_sat_rel_earth = r_sat_rel_moon' + r_moon;
    v_sat_rel_earth = v_sat_rel_moon' + v_moon;

    % Convert position and velocity vectors back to orbital elements
    [a, ecc, incl, RAAN, argp, nu] = ijk2keplerian(r_sat_rel_earth*1e3, v_sat_rel_earth*1e3);   % m,m/s

end