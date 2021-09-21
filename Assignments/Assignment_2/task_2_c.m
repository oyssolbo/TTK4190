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

% Rewriting into Evans form
sys = tf([a_phi_2], [1 a_phi_1+a_phi_2*K_d a_phi_2*K_p 0]);

K_min = -3.2;
K_max = 0;
dk = 0.1;
K = [K_min:dk:K_max];

%h = rlocus(sys,K);
%p = getoptions(h);
%p.Title.String = 'Locus plot for K = [-3.2, 0]';
%setoptions(h, p);

rlocus_handle = rlocusplot(sys, K);
rlocus_opt = getoptions(rlocus_handle);
rlocus_opt.Title.String = 'Locus plot for K = [-3.2, 0]';
setoptions(rlocus_handle, rlocus_opt);