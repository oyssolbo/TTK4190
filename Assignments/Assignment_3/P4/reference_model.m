function [r_d_out, a_d_out, psi_d_out] ...
    = reference_model(r_d, a_d, psi_d, psi_ref, h)

    w_ref = 0.03;
    zeta_ref = 1;

    r_max = 1 * pi/180;
    a_max = 0.5 * pi/180;

    sat_rd = sat_value(r_d, r_max);
    sat_ad = sat_value(a_d, a_max);
    
    psi_d_dot = sat_rd;
    r_d_dot = sat_ad;
    a_d_dot = -(2*zeta_ref + 1)*w_ref*sat_ad ...
              -(2*zeta_ref + 1)*w_ref^2*sat_rd ...
              + w_ref^3*(psi_ref - psi_d);
    
    psi_d_out = psi_d + h*psi_d_dot;
    r_d_out = r_d + h*r_d_dot;
    a_d_out = a_d + h*a_d_dot;

end

function sat_val = sat_value(val, abs_lim)
    sat_val = val;
    if abs(val) >= abs_lim
       sat_val = sign(val)*abs_lim; 
    end
end
