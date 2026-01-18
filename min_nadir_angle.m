function theta_deg = min_nadir_angle(pos_sat, pos_body, r_body)

    dist_to_center = norm(pos_body-pos_sat);
    theta = asin(r_body/dist_to_center);
    theta_deg = rad2deg(theta);

end