%% Description
% M-script for numerically solving the differential equations for a course
% controller. The course controller consists of a controller in cascade
% with a roll controller
%
% System equations
%   p_dot = -a_phi_1*p + a_phi_2*delta_a
%
%   chi_dot = g/V_g*(phi + d)
%
%   delta_a_dot = -T_a*delta_a + T_a*delta_a^c
%
%   phi_ref = e_chi*(Ki/s + Kp)
%
%   delta_a_ref = e_phi*Kp + k_d*p

%%
clear;
clc;
close all;

%% Simulation time
h = 0.001;                     % sample time (s)
N  = 5000;                   % number of samples. Should be adjusted
t = 1:N+1;

%% Constants
deg2rad = pi/180;   
rad2deg = 180/pi;

g = 9.81;

%% System
% x = [chi, beta, phi, p, r, delta_a]^T

% Countinous models
A_c = [-0.322   0.052   0.028  -1.12    0.002;
        0       0       1      -0.001   0;
      -10.6     0      -2.87    0.46   -0.65;
        6.87    0      -0.04   -0.32   -0.02;
        0       0       0       0      -7.5];
   
B_c = [0;
       0;
       0;
       0;
       7.5];
 
C_c = [1 0 0 0 0;
       0 1 0 0 0;
       0 0 1 0 0;
       0 0 0 1 0];
   
D_c = zeros(4,1);

%% Model parameters
V_g = 540/3.6; % SOG (assuming V_w = 0) in m/s

a_phi_1 = 2.87;
a_phi_2 = -0.65;

T_a = 7.5;
W_chi = 10; % Should be between 5 and 10

delta_a_max = 30;   % deg
e_phi_max = 15;     % deg
zeta_phi = sqrt(2)/2; % ~0.707

zeta_chi = 1; % > 0

% Exploiting eq. 6.7 to 6.9 in B&W
k_p_phi = delta_a_max/e_phi_max * sign(a_phi_2);
w_phi_n = sqrt(abs(a_phi_2) * delta_a_max/e_phi_max);
k_d_phi = (2*zeta_phi*w_phi_n - a_phi_1)/a_phi_2;

w_n_chi = 1/W_chi * w_phi_n;

% Exploiting 6.12 & 6.12 in B&W 
k_p_chi = 2*zeta_chi*w_n_chi*V_g/g;
k_i_chi = w_n_chi^2*V_g/g;

%% Disturbance
d_chi = 1.5*deg2rad;

%% Initial and reference values
chi_0 = ssa(350)*deg2rad;
beta = 0;
phi = ssa(-10)*deg2rad;
p = 0; 
r = 0;
delta_a = 0;

x = [beta; phi; p; r; delta_a]; 

chi_ref = zeros(N+1, 1);
chi_ref(N/5:N/2) = 15*deg2rad;
chi_ref(N/2+1:floor(2*N/3)) = 0;
chi_ref(floor(2*N/3)+1:N+1) = -15*deg2rad;
%chi_ref(N/2+1:N+1) = -15*deg2rad; % Testing purposes only

% Creating a time-series
chi_ref_ts = timeseries(chi_ref, t);

% Trying to circumwent problems with timeseries
chi_ref_step = [zeros(1,N/5) 15*ones(1,N/2-N/5)*deg2rad zeros(1,floor(2/3*N)-N/2) -15*ones(1,N-floor(2/3*N))*deg2rad];
chi_ref_step = [zeros(1,10), 15*deg2rad*ones(1,10), zeros(1,10), -15*deg2rad*ones(1,10)];
e_chi_int = 0; % The integrator for chi
%% Memory allocation
table = zeros(N+1,5); % chi, phi, delta_a, p , phi_ref   

sim("simulink_2f.slx");
load("simulation_results.mat");

%% Simulation loop
%{
for i = 1:N+1,
    dt = (i-1)*h; 
    
    % Transforming to 0-360 deg
    chi = ssa(chi);
    
    % Limit the gain in delta_a
    if(abs(x(5))*rad2deg >= delta_a_max)
       x(5) = sign(x(5))*delta_a_max*deg2rad; 
    end
    
    % Calculate system states
    e_chi = ssa(chi_ref(i) - chi);
    phi_ref = k_p_chi*e_chi + k_i_chi*e_chi_int;
    e_phi = ssa(phi_ref - x(2));
    delta_a_ref = e_phi*k_p_phi + k_d_phi*x(3);
    
    % Calculate the differential equations
    chi_dot = g/V_g*tan(x(2))*cos(chi-x(2)) + d_chi;
    e_chi_int_dot = e_chi;
    x_dot = A_c*x + B_c*delta_a_ref;
    
    % Store values of step i
    table(i,:) = [chi x(2) x(5) x(3) phi_ref]; 
    
    % Calculate values for i+1
    x = x + h*x_dot;
    chi = chi + h*chi_dot;
    e_chi_int = e_chi_int + h*e_chi_int_dot;
end 
%}

%% Run simulink simulation

%% PLOT FIGURES
t         = ans(1,:);
chi       = rad2deg*ans(2,:);  
phi       = rad2deg*ans(3,:); 
delta_a   = rad2deg*ans(4,:);
p         = rad2deg*ans(5,:);
phi_ref   = rad2deg*ans(6,:);
chi_ref   = ans(7,:);

figure (1); clf;
hold on;
plot(t, phi);
plot(t, phi_ref);
hold off;
grid on;
legend('\phi', '\phi_{ref}');
title('Roll');
xlabel('Time [s]'); 
ylabel('Angle [deg]');

figure (2); clf;
hold on;
plot(t, chi);
plot(t, chi_ref*rad2deg);
hold off;
grid on;
legend('\chi', '\chi_{ref}');
title('Course');
xlabel('Time [s]'); 
ylabel('Angle [deg]');

figure (3); clf;
hold on;
plot(t, delta_a);
hold off;
grid on;
legend('\delta_a');
title('Rudder input');
xlabel('Time [s]'); 
ylabel('Angle [deg]');

figure (4); clf;
hold on;
plot(t, p);
hold off;
grid on;
legend('p');
title('Roll rate');
xlabel('Time [s]'); 
ylabel('Angular velocity [deg/s]');
