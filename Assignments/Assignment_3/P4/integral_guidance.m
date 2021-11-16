function [chi_d, dy_int] = integral_guidance(pos, wp_ref, wp_t, delta, kappa, y_int)
    % pos = (x, y) position of ownship in NED
    % wp_ref = reference waypoint, waypoint 1
    % wp_t = target waypoint, waypoint 2
    % lookahead = tuning parameter
    
    pi_p = atan2(wp_t(2) - wp_ref(2), wp_t(1) - wp_ref(1));
    [~, ~, y_e] = crosstrack(wp_t(1), wp_t(2), ...
                             wp_ref(1), wp_ref(2), ...
                             pos(1), pos(2));
    Kp = 1 / delta;
    Ki = kappa * Kp;         
                         
    chi_d = pi_p - atan(Kp*y_e + Ki*y_int);
    dy_int = delta*y_e / (delta^2 + (y_e + kappa*y_int)^2);
end