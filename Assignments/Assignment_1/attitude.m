% M-script for numerical integration of the attitude dynamics of a rigid 
% body represented by unit quaternions. The MSS m-files must be on your
% Matlab path in order to run the script.
%
% System:                      .
%                              q = T(q)w
%                              .
%                            I w - S(Iw)w = tau
% Control law:
%                            tau = constant
% 
% Definitions:             
%                            I = inertia matrix (3x3)
%                            S(w) = skew-symmetric matrix (3x3)
%                            T(q) = transformation matrix (4x3)
%                            tau = control input (3x1)
%                            w = angular velocity vector (3x1)
%                            q = unit quaternion vector (4x1)
%
% Author:                   2018-08-15 Thor I. Fossen and Håkon H. Helgesen

%%
clear;
clc;
close all;

%% USER INPUTS
h = 0.1;                     % sample time (s)
N = 1000;                    % number of samples. Should be adjusted

% Model parameters
m = 180;
r = 2;
I = m*r^2*eye(3);            % inertia matrix
I_inv = inv(I);

% Constants
deg2rad = pi/180;   
rad2deg = 180/pi;

phi = -5*deg2rad;            % initial Euler angles
theta = 10*deg2rad;
psi = -20*deg2rad;

q = euler2q(phi,theta,psi);   % transform initial Euler angles to q

w = [0 0 0]';                 % initial angular rates

table = zeros(N+1,14);        % memory allocation
tracking_error = zeros(N+1, 3);
ref_signal = zeros(N+1, 3);

% Control law
k_p = 20;
k_d = 400;

K_d = k_d * eye(3);

% Desired system state
t = 1:N + 1;
phi_t = 0 * t * deg2rad;
theta_t = 15 * cos(0.1 * t) * deg2rad;
psi_t = 10 * sin(0.05 * t) * deg2rad;

euler_d = [phi_t; theta_t; psi_t];
euler_dot_d = [0 * t * deg2rad;
               -1.5 * sin(0.1 * t) * deg2rad;
               0.5 * cos(0.05 * t) * deg2rad];

eta_t = sqrt(1 - phi_t.^2 - theta_t.^2 - psi_t.^2);

%q_d = cell2mat(arrayfun(@(phi, theta, psi) euler2q(phi, theta, psi), ...
%    phi_t, theta_t, psi_t, 'UniformOutput', false))';
q_d = [eta_t; phi_t; theta_t; psi_t]';
q_d_conj = quatconj(q_d);

w_d = [0, 0, 0]';

%% FOR-END LOOP
for i = 1:N+1,
    t = (i-1)*h;                  % time
   
    [phi,theta,psi] = q2euler(q); % transform q to Euler angles
    [J,J1,J2] = quatern(q);       % kinematic transformation matrices
   
    q_dot = J2*w;                        % quaternion kinematics
   
    e = q(2:4);
    
    q_tilde = quatmultiply(q_d_conj(i,:), q');
    e_tilde = q_tilde(2:4);
    tracking_error(i,:) = e_tilde;
    
    T_euler_inv = [1, 0, -sin(theta);
              0, cos(phi), cos(theta) * sin(phi);
              0, -sin(phi), cos(theta) * cos(phi)];
    
    % Control law for problem 1.2
    tau = control_law(e, w, K_d, k_p);
   
    % Control law for problem 1.5 and 1.6 (quaternion)
    tau = control_law(e_tilde', w, K_d, k_p);
   
    % Control law for problem 1.6 (angular velocity)
    w_d = T_euler_inv * euler_dot_d(:,i);
    w_tilde = w - w_d;
    tau = control_law(e_tilde', w_tilde, K_d, k_p);
   
    %q_dot = J2*w;                        % quaternion kinematics
    w_dot = I\(Smtrx(I*w)*w + tau);  % rigid-body kinetics
   
    table(i,:) = [t q' phi theta psi w' tau'];  % store data in table
   
    q = q + h*q_dot;	             % Euler integration
    w = w + h*w_dot;
   
    q  = q/norm(q);               % unit quaternion normalization
end 

%% PLOT FIGURES
t       = table(:,1);  
q       = table(:,2:5); 
phi     = rad2deg*table(:,6);
theta   = rad2deg*table(:,7);
psi     = rad2deg*table(:,8);
w       = rad2deg*table(:,9:11);  
tau     = table(:,12:14);

e_tilde = tracking_error(:,1:3) .* rad2deg;

figure(1);
brg_plot(t, [phi theta psi]);
grid on;
legend('\phi', '\theta', '\psi');
title('Euler angles');
xlabel('time [s]'); 
ylabel('angle [deg]');

figure(2);
brg_plot(t, w);
grid on;
legend('p', 'q', 'r');
title('Angular velocities');
xlabel('time [s]'); 
ylabel('angular rate [deg/s]');

figure(3);
brg_plot(t, tau);
grid on;
legend('x', 'y', 'z');
title('Control input');
xlabel('time [s]'); 
ylabel('input [Nm]');

figure(4);
brg_plot(t, e_tilde);
grid on;
legend('e_1', 'e_2', 'e_3');
title('Tracking error - epsilon');
xlabel('time [s]'); 
ylabel('angle [deg]');

figure(5);
brg_plot(t, ref_signal);
grid on;
legend('e_1', 'e_2', 'e_3');
title('Reference signal');
xlabel('time [s]'); 
ylabel('angle [deg]');

function tau = control_law(e, w, K_d, k_p)
    tau = - K_d * w - k_p * e;
end

function brg_plot(t, v)
    hold on;
    plot(t, v(:,1), 'b');
    plot(t, v(:,2), 'r');
    plot(t, v(:,3), 'g');
    hold off;
end