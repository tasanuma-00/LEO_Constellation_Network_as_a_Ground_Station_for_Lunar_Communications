function theta_deg = angle_between_vectors(A, B)
    % Calculate the dot product and magnitudes of the vectors
    dotProd = dot(A, B);
    magA = norm(A);
    magB = norm(B);

    % Calculate the angle in radians
    theta = acos(dotProd / (magA * magB));

    % Convert to degrees
    theta_deg = rad2deg(theta);
end