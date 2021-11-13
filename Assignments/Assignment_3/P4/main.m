% Project in TTK4190 Guidance, Navigation and Control of Vehicles 
%
% Author:           Solbø, Ø. & Strøm, C.
% Study program:    MTTK

clear; 
close all;
%% Input
h  = 0.1;               % Sampling time [s]
Ns = 100000;             % Num samples
U_ref = 9;            % Desired surge speed (m/s)

%% Initial states
eta = [0 0 deg2rad(-110)]';
nu  = [0.1 0 0]';
delta = 0;
n = 0;
Qm = 0;
x = [nu' eta' delta n Qm]';

%% Ship coefficients
rudder_max  = 40 * pi/180;        % Max rudder angle      (rad)
d_rudder_max = 5  * pi/180;       % Max rudder derivative (rad/s)

Loa = 161;

% First order Nomoto
K = 0.0075;
T = 169.5493;

%% Guidance
load WP
lookahead = 15*Loa;
ROA = 5*Loa;
wp_n = 1;

%% External forces
Vc = 0;
beta_vc = deg2rad(45);

% Wind-coefficients
Vw = 0;
beta_vw = deg2rad(135);
rho_a = 1.247;
cy = 0.95;
cn = 0.15;
A_Lw = 10*Loa;

% Initial wind-force/moment
Ywind = 0;
Nwind = 0;

%% Controller
zeta_n = 1;
w_b = 0.06; 
w_n = 1/(sqrt(1-2*zeta_n^2 + sqrt(4*zeta_n^4 - 4*zeta_n^2 + 2))) * w_b;

% Using eq. 15.95, ex. 15.7 and algorithm 15.1
m = T/K;
d = 1/K;

Kp = m*w_n^2;
Kd = 2*zeta_n*w_n*m - d;
Ki = w_n/10*Kp;

%% Reference model
psi_d = 0;
r_d = 0;
a_d = 0;

e_psi_int = 0;
delta_c = 0;
n_c = 9;

%% Simulation
simdata = zeros(Ns+1,18);       % Table of simulation data

for i=1:Ns+1
    %% Time
    t = (i-1) * h;              % Time (s)
    
    %% Currents
    nu_c = [Vc*cos(beta_vc - x(6)), Vc*sin(beta_vc - x(6)), 0 ]';
    
    %% Wind
    if t >= 200
        u_rw = x(1) - Vw*cos(beta_vw - x(6));
        v_rw = x(2) - Vw*sin(beta_vw - x(6));
        V_rw = norm([u_rw, v_rw], 2);
        
        gamma_rw = -atan2(u_rw, v_rw);
        dp = 1/2*rho_a*V_rw^2;
        
        Ywind = dp*cy*sin(gamma_rw)*A_Lw;
        Nwind = dp*cn*sin(2*gamma_rw)*A_Lw*Loa;
    end
    tau_wind = [0 Ywind Nwind]';
    
    %% Guidance and reference model
    wp_ref = WP(:, wp_n);
    wp_t = WP(:, wp_n+1);
    
    % Manage waypoint switching
    if norm(wp_t - x(4:5), 2) <= ROA
        wp_n = wp_n + 1;
        if wp_n == width(WP)
            break % End simulation if we reach the final waypoint
        end
    end
    
    
  
    psi_ref = guidance(x(4:5), wp_ref, wp_t, lookahead);
    psi_ref = wrapTo2Pi(psi_ref);
    
    [r_d, a_d, psi_d] = reference_model(r_d, a_d, psi_d, psi_ref, h);
    
    %% Control law
    % Heading
    e_psi = ssa(-psi_d + x(6));
    e_r = ssa(-r_d + x(3));
    
    % With integral windup
    e_psi_int_dot = anti_windup(e_psi, delta_c, delta_c - x(7), rudder_max, d_rudder_max);
    e_psi_int = e_psi_int + h*e_psi_int_dot;
    
    delta_c = -Kp*e_psi - Kd*e_r - Ki*e_psi_int;
    %% Ship dynamics
    u = [delta_c n_c]';
    [xdot,u] = ship(x,u,nu_c,tau_wind);
    
    %% Sideslip and crab
    crab = atan2(x(2), x(1));
    
    R = Rzyx(0, 0, x(6));
    U_ned = R'*[x(1), x(2), 0]';
    U_r = U_ned - nu_c';
    v_r = U_r(2);
    
    sideslip = asin(v_r/norm(U_r, 2));
    
    %% Course
    chi = wrapTo2Pi(x(6) + crab);
    chi_ref = wrapTo2Pi(psi_ref);
    
    %% Store simulation data 
    simdata(i,:) = [t x(1:3)' x(4:6)' x(7) x(8) u(1) u(2) U_ref psi_d r_d crab sideslip, chi, chi_ref];     
 
    %% Euler integration
    x = euler2(xdot,x,h);  
end

simdata( all(~simdata,2), : ) = []; % Cut trailing zeros

%% Plotting
t       = simdata(:,1);                 % s
u       = simdata(:,2);                 % m/s
v       = simdata(:,3);                 % m/s
r       = (180/pi) * simdata(:,4);      % deg/s
x       = simdata(:,5);                 % m
y       = simdata(:,6);                 % m
psi     = (180/pi) * wrapTo2Pi(simdata(:,7));      % deg
delta   = (180/pi) * simdata(:,8);      % deg
n       = 60 * simdata(:,9);            % rpm
delta_c = (180/pi) * simdata(:,10);     % deg
n_c     = 60 * simdata(:,11);           % rpm
u_d     = simdata(:,12);                % m/s
psi_d   = (180/pi) * simdata(:,13);     % deg
r_d     = (180/pi) * simdata(:,14);     % deg/s
crab    = (180/pi) * simdata(:,15);     % deg
sideslip= (180/pi) * simdata(:,16);     % deg
chi     = (180/pi) * simdata(:,17);     % deg
chi_d = (180/pi) * simdata(:,18);     % deg

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

figure(5)
plot(t,chi,t,chi_d,'linewidth',2);
title('Actual and desired course (deg)');
legend({'\chi', '\chi_d'})
xlabel('Time (s)');
ylabel('Course (deg)');

pathplotter(x,y);

%% Functions
function x_dot_aw = anti_windup(x_dot, param_0_check, param_1_check , abs_0, abs_1)
    x_dot_aw = x_dot;
    if abs(param_0_check) >= abs_0 || abs(param_1_check) >= abs_1
        x_dot_aw = 0; 
    end
end