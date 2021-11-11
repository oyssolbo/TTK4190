function [chi_d, y_int_dot] = guidance(state_struct, current_waypoint, next_waypoint, ILOS_params)
    % Important note regarding adpative lookahead distance
    % A larger lookadhead distance is better when the system is closer to
    % the desired waypoint, as it prevents overshoot better
    % A smaller lookahead distance is preferable when the ships is further
    % away, as it creates a more aggressive controller
    
    % Assumptions:
    % - The values for position and waypoints are in NED, such that this
    % function must be modified if used over larger areas
    % - The function should be invoked for every iteration, or often
    % enough, such that will generate a new desired course whenever one
    % detects that the old desired course no longer is satisfying. Since a
    % ship is a slow system, the function could be invoked for every tenth
    % iteration or something
    % - The funstion takes in the lookahead distance, such that this could
    % be chosen spesifically between ships of different types. The same
    % reason is also chosen for kappa
    
    % Extract values
    lookahead = ILOS_params.lookahead;
    kappa = ILOS_params.kappa;
    
    pos = state_struct.pos;
    y_int = state_struct.y_int;
    
    % Calculate pi_p
    diff_waypoints = next_waypoint - current_waypoint;
    pi_p = atan2(diff_waypoints(2), diff_waypoints(1));
    
    % Using eq. 12.55 to calculate the cross-track error
    R = [cos(pi_p), -sin(pi_p);
         sin(pi_p),  cos(pi_p)];
    ep = R' * (pos - current_waypoint);
    y_ep = ep(2);
    
    % Calculate controller gains according to 12.79 and 12.109
    Kp = 1/lookahead;
    Ki = kappa*Kp;
    
    % Calculate desired course using eq. 12.78
    chi_d = pi_p - atan(Kp * y_ep + Ki*y_int);
    
    % Calculate y_int_dot using eq. 12.109
    % Should ideally be calculated outside of this function, but the
    % consistency of this code is so terrible nonetheless...
    y_int_dot = lookahead*y_ep / (lookahead^2 + (y_ep + kappa*y_int)^2);
end