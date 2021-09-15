% Assignment 2 
% Task 2c)

a_phi_1 = 2.87;
a_phi_2 = -0.65;

zeta = 0.707;
e_max = 15;
delta_max = 30;

% According to 6.7 to 6.9 in B&W
K_p = delta_max / e_max * sign(a_phi_2);
w_n = sqrt(abs(a_phi_2) * delta_max / e_max);
K_d = (2 * zeta * w_n - a_phi_1) / a_phi_2;

% num = a_phi_2*K_p*s + K_i;
% den = s^3 + (a_phi_1+a_phi_2*K_d)*s^2 + a_phi_2*K_p*s + K_i*a_phi_2;  

%sys = num / den;
%sys = tf([a_phi_2*K_p K_i],[1 a_phi_1+a_phi_2*K_d a_phi_2*K_p K_i*a_phi_2])

% Rewriting into Evans form
sys = tf([a_phi_2], [1 a_phi_1+a_phi_2*K_d a_phi_2*K_p 0])

rlocus(sys, [-3.2:0.1:0]);

% Can see that approximating Ki in range [-3.2, 0] is required for
% the system to be stable AKA have poles in the left half plane


