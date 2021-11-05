% Project in TTK4190 Guidance, Navigation and Control of Vehicles 
%
% Author:           Solbø, Ø. & Strøm, C.
% Study program:    MTTK

clear all; close all;
%% Input
h  = 0.1;               % Sampling time [s]
Ns = 10000;             % Num samples

psi_ref_init = 10 * pi/180 * ones(1, Ns/4);  % Desired yaw angle (rad)
psi_ref_end = -20 * pi/180 * ones(1, (Ns - Ns/4) + 1);
psi_ref = [psi_ref_init, psi_ref_end];
%psi_ref = 10 * pi/180;
U_ref   = 9;            % Desired surge speed (m/s)

%% Initial states
eta = [0 0 0]';
nu  = [0.1 0 0]';
delta = 0;
n = 0;
x = [nu' eta' delta n]';

%% Ship coefficients
% Rudder limitations
rudder_max  = 40 * pi/180;        % Max rudder angle      (rad)
d_rudder_max = 5  * pi/180;       % Max rudder derivative (rad/s)

Loa = 161;
m = 17.0677e6;

% First order Nomoto
K = 0.0075;
T = 169.5493;

% Thrust-calculation
J_a = 0;    % Open-water advance coefficient
PD = 1.5;   % Pitch / diameter ratio
AEAO = 0.65;% Blade area ratio
z = 4;      % Num propeller-blades
[K_T, K_Q] = wageningen(J_a, PD, AEAO, z);

%% Propultion system
I_m = 1e5;
K_m = 0.6;
T_m = 10;
tau = 0;
rho = 1025;
epsilon = 1e-14;

%J_a = @(n, u) u_a / (n * D + epsilon);  % Adding epsilong to prevent division by zero
%T = @(n, u_a) rho*D^4*K_T*J_a(n, u)*abs(n)*n;
%Q = @(n, u_a) rho*D^5*K_Q*J_a(n, u)*abs(n)*n;

%% External forces
Vc = 1;
beta_vc = deg2rad(45);

% Wind-coefficients
Vw = 10;
beta_vw = deg2rad(135);
rho_a = 1.247;
cy = 0.95;
cn = 0.15;
A_Lw = 10*Loa;

%% Controller
zeta_PID = 1;
w_b = 0.06; 

% Calculating controller-gains using p. 537
% The controller becomes unstable due to having too large controller
% valus due the large mass m
w_PID = 1/(sqrt(1-2*zeta_PID^2 + sqrt(4*zeta_PID^4 - 4*zeta_PID^2 + 2))) * w_b;

Kp = m*w_PID^2;
Kd = 2*zeta_PID*w_PID*m;
Ki = w_PID/10*Kp;

% Using eq. 15.69 - 15.71 to calculate the controller gains. Gives a much
% better response compared to preiously
e_psi_max = 1 * pi/180;
delta_r_max = 40 * pi/180;

w_psi = sqrt(K/T * delta_r_max/e_psi_max);

Kp = delta_r_max / e_psi_max;
Kd = (2*zeta_PID*w_psi*T - 1)/K;
Ki = w_psi/10 * Kp;

%% Reference model
w_ref = 0.03;
zeta_ref = 1;

psi_d = 0;
r_d = 0;
a_d = 0;

% But shouldn't these also be used to represent physical limitations
% that is present in the ship/model?
r_max = 1 * pi/180;
a_max = 0.5 * pi/180;

e_psi_int = 0;
delta_c = 0;

%% Simulation
simdata = zeros(Ns+1,16);       % Table of simulation data
psi_d_arr = zeros(1, Ns+1);

for i=1:Ns+1
    %% Time
    t = (i-1) * h;              % Time (s)
    
    %% Current disturbance
    uc = Vc*cos(beta_vc - x(6));
    vc = Vc*sin(beta_vc - x(6));
    nu_c = [ uc vc 0 ]';
    
    %% Wind disturbance
    Ywind = 0;
    Nwind = 0;
    if t >= 200000,
        u_rw = x(1) - Vw*cos(beta_vw - x(6));
        v_rw = x(2) - Vw*sin(beta_vw - x(6));
        V_rw = norm([u_rw, v_rw], 2);
        gamma_rw = -atan2(u_rw, v_rw);
        
        dynamic_pressure = 1/2*rho_a*V_rw^2;
        Ywind = dynamic_pressure*cy*sin(gamma_rw)*A_Lw;
        Nwind = dynamic_pressure*cn*sin(2*gamma_rw)*A_Lw*Loa;
    end
    tau_wind = [0 Ywind Nwind]';
    
    u_d = U_ref;
    %delta_c = 7.5;
    %n_c = 10;
    
    %% Reference model
    sat_rd = sat_value(r_d, r_max);
    sat_ad = sat_value(a_d, a_max);
    psi_d_dot = sat_rd;
    r_d_dot = sat_ad;
    a_d_dot = -(2*zeta_ref + 1)*w_ref*sat_ad - (2*zeta_ref + 1)*w_ref^2*sat_rd + w_ref^3*(psi_ref(i) - psi_d);   % -(2*zeta_ref + 1)*w_ref*sat(a_d) - (2*zeta_ref + 1)*w_ref^2*sat(r_d) + w_ref^3*(psi_ref - psi_d);
    
    psi_d = psi_d + h*psi_d_dot;
    r_d = r_d + h*r_d_dot;
    a_d = a_d + h*a_d_dot;
    
    psi_d_arr(i) = psi_d;
    
    %% Control law
    % Heading
    e_psi = ssa(-psi_d + x(6));
    e_r = ssa(-r_d + x(3));
    
    % With integral windup
    e_psi_int_dot = anti_windup(e_psi, delta_c, delta_c - x(7), rudder_max, d_rudder_max);
    e_psi_int = e_psi_int + h*e_psi_int_dot;
    
    delta_c = -Kp*e_psi - Kd*e_r - Ki*e_psi_int;
    
    % Propeller speed
    n_c = 10;                   % Propeller speed (rps) (But is it really in rps or in rpm as indicated in ship.m)
    
    %% Ship dynamics
    u = [delta_c n_c]';
    [xdot,u] = ship(x,u,nu_c,tau_wind, U_ref);
    
    %% Sideslip and crab
    % Using equation 10.138, it is given that the body-velocities are
    % calculated in the ship dynamics
    % Calculating crab and sideslip based on this
    crab = atan2(x(2), x(1));
    
    R = Rzyx(0, 0, x(6));
    U_ned = R'*[x(1), x(2), 0]';
    U_r = U_ned - [uc, vc, 0]';
    u_r = U_r(1);
    
    sideslip = asin(u_r/norm(U_r, 2));
    
    %% Store simulation data 
    simdata(i,:) = [t x(1:3)' x(4:6)' x(7) x(8) u(1) u(2) u_d psi_d r_d crab sideslip];     
 
    %% Euler integration
    x = euler2(xdot,x,h);    
end

%% Plotting
t       = simdata(:,1);                 % s
u       = simdata(:,2);                 % m/s
v       = simdata(:,3);                 % m/s
r       = (180/pi) * simdata(:,4);      % deg/s
x       = simdata(:,5);                 % m
y       = simdata(:,6);                 % m
psi     = (180/pi) * simdata(:,7);      % deg
delta   = (180/pi) * simdata(:,8);      % deg
n       = 60 * simdata(:,9);            % rpm
delta_c = (180/pi) * simdata(:,10);     % deg
n_c     = 60 * simdata(:,11);           % rpm
u_d     = simdata(:,12);                % m/s
psi_d   = (180/pi) * simdata(:,13);     % deg
r_d     = (180/pi) * simdata(:,14);     % deg/s
crab    = (180/pi) * simdata(:,15);     % deg
sideslip= (180/pi) * simdata(:,16);     % deg

figure(1)
figure(gcf)
subplot(311)
plot(y,x,'linewidth',2); axis('equal')
title('North-East positions (m)'); 
xlabel('(m)'); ylabel('(m)'); 
subplot(312)
plot(t,psi,t,psi_d,'linewidth',2);
title('Actual and desired yaw angles (deg)'); 
xlabel('time (s)');
subplot(313)
plot(t,r,t,r_d,'linewidth',2);
title('Actual and desired yaw rates (deg/s)'); 
xlabel('time (s)');

figure(2)
figure(gcf)
subplot(311)
plot(t,u,t,u_d,'linewidth',2);
title('Actual and desired surge velocities (m/s)'); 
xlabel('time (s)');
subplot(312)
plot(t,n,t,n_c,'linewidth',2);
title('Actual and commanded propeller speed (rpm)'); 
xlabel('time (s)');
subplot(313)
plot(t,delta,t,delta_c,'linewidth',2);
title('Actual and commanded rudder angles (deg)'); 
xlabel('time (s)');

figure(3)
figure(gcf)
subplot(211)
plot(t,crab,'linewidth',2);
title('Crab (deg)'); xlabel('time (s)');
subplot(212)
plot(t,sideslip,'linewidth',2);
title('Sideslip (deg)'); xlabel('time (s)');

figure(4)
plot(t,u,t,u_d,'linewidth',2);
title('Actual and desired surge velocities (m/s)');
legend({'u', 'u_d'})
xlabel('Time (s)');
ylabel('Surge (m/s)');

%% Functions
function sat_val = sat_value(val, abs_lim)
    sat_val = val;
    if abs(val) >= abs_lim
       sat_val = sign(val)*abs_lim; 
    end
end

function x_dot_aw = anti_windup(x_dot, param_0_check, param_1_check , abs_0, abs_1)
    x_dot_aw = x_dot;
    if abs(param_0_check) >= abs_0 || abs(param_1_check) >= abs_1
        x_dot_aw = 0; 
    end
end