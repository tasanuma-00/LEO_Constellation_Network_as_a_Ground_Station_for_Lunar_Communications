function is_accessible = is_leo_sat_accessible(leo_sat_pos,target_pos,earth_pos,cone_angle_threshold)

% Moon radius in km
R_earth = 6371;       

% set accessible to false by default
is_accessible = false;

% initialize vector
leo_sat_nadir_vector = earth_pos - leo_sat_pos;
leo_sat_zenith_vector = leo_sat_pos - earth_pos;
leo_sat_to_target_vector = target_pos - leo_sat_pos;

% calculate angles for checking blockage with Earth
angle_separation = angle_between_vectors(leo_sat_nadir_vector,leo_sat_to_target_vector);
earth_blockage_angle = min_nadir_angle(leo_sat_pos, earth_pos, R_earth*1e3);

% calculate angles for checking cone angle within bound
leo_sat_cone_angle = angle_between_vectors(leo_sat_zenith_vector,leo_sat_to_target_vector);

% check if both conditions met
if angle_separation > earth_blockage_angle && leo_sat_cone_angle < cone_angle_threshold
    is_accessible = true;
end

end