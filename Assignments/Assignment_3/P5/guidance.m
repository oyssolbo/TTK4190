function chi_d = guidance(pos, wp_ref, wp_t, lookahead)
    % pos = (x, y) position of ownship in NED
    % wp_ref = reference waypoint, waypoint 1
    % wp_t = target waypoint, waypoint 2
    % lookahead = tuning parameter
    
    pi_p = atan2(wp_t(2) - wp_ref(2), wp_t(1) - wp_ref(1));
    [~, ~, y_e] = crosstrack(wp_t(1), wp_t(2), ...
                             wp_ref(1), wp_ref(2), ...
                             pos(1), pos(2));
    chi_d = pi_p - atan(y_e / lookahead);
end