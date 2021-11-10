function chi_d = guidance(ship_pos, current_waypoint, next_waypoint, lookahead)
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
    % be chosen spesifically between ships of diferent types
    
    % Calculate pi_p
    diff_waypoints = next_waypoint - current_waypoint;
    pi_p = atan2(diff_waypoints(2), diff_waypoints(1));
    
    % Using eq. 12.55 to calculate the cross-track error
    R = [cos(pi_p), -sin(pi_p);
         sin(pi_p),  cos(pi_p)];
    ep = R' * (ship_pos - current_waypoint);
    y_ep = ep(2);
   
    % Calculate controller gains according to 12.79
    Kp = 1/lookahead;
    
    % Calculate desired course uisng eq. 12.78
    chi_d = pi_p - atan(Kp * y_ep);
end