function is_accessible = is_lunar_sat_accessible(sat_pos,target_pos,moon_pos,cone_angle_threshold)

% Moon radius in km
R_moon = 1737.4;       

% set accessible to false by default
is_accessible = false;

% initialize vector
lunar_sat_nadir_vector = moon_pos - sat_pos;
lunar_sat_to_target_vector = target_pos - sat_pos;
lunar_sat_zenith_vector = sat_pos - moon_pos;

% calculate angles for checking blockage with moon
angle_separation = angle_between_vectors(lunar_sat_nadir_vector,lunar_sat_to_target_vector);
moon_blockage_angle = min_nadir_angle(sat_pos, moon_pos, R_moon*1e3);

% calculate angles for checking cone angle within bound
lunar_sat_cone_angle = angle_between_vectors(lunar_sat_zenith_vector,lunar_sat_to_target_vector);

% check if both conditions met
if angle_separation > moon_blockage_angle && lunar_sat_cone_angle < cone_angle_threshold
    is_accessible = true;
end

end