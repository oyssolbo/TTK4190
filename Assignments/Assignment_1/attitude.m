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
%clear;
%clc;
%close all;

%% USER INPUTS
h = 0.1;                     % sample time (s)
N = 5000;                    % number of samples. Should be adjusted

% Model parameters
m = 180;
r = 2;
I = m*r^2*eye(3);            % inertia matrix
I_inv = inv(I);

% Control law
k_p = 20;
k_d = 400;
K_d = k_d * eye(3);

% Constants
deg2rad = pi/180;   
rad2deg = 180/pi;

phi = -5*deg2rad;            % initial Euler angles
theta = 10*deg2rad;
psi = -20*deg2rad;

q = euler2q(phi,theta,psi);   % transform initial Euler angles to q

w = [0 0 0]';                 % initial angular rates
w_d = [0, 0, 0]';

table = zeros(N+1,20);        % memory allocation
task = "1.6";
%% FOR-END LOOP
for i = 1:N+1
    t = (i-1)*h;                  % time
   
    [phi,theta,psi] = q2euler(q); % transform q to Euler angles
    [J,J1,J2] = quatern(q);       % kinematic transformation matrices
    
    phi_d   = 0;
    theta_d = 15 * cos(0.1  * t) * deg2rad;
    psi_d   = 10 * sin(0.05 * t) * deg2rad;
    
    dTheta_d = [0;
               -1.5 * sin(0.1  * t) * deg2rad;
                0.5 * cos(0.05 * t) * deg2rad];
        
    q_d     = euler2q(phi_d, theta_d, psi_d);
    q_d_bar = quatconj(q_d');
 
    q_tilde = quatmultiply(q_d_bar, q');
    e_tilde = q_tilde(2:4);
   
    % Control law
    switch task
        case "1.2"
            tau = control_law(e, w, K_d, k_p);
        case "1.5"
            tau = control_law(e_tilde', w, K_d, k_p);
        case "1.6"
            w_tilde = w - Tzyx(phi_d, theta_d) \ dTheta_d;
            tau = control_law(e_tilde', w_tilde, K_d, k_p);
    end

    table(i,:) = [t q' phi theta psi w' tau' phi_d theta_d psi_d e_tilde];  % store data in table
    
    q_dot = J2*w;                    % quaternion kinematics
    w_dot = I\(Smtrx(I*w)*w + tau);  % rigid-body kinetics
    
    q = q + h*q_dot;	             % Euler integration
    w = w + h*w_dot;
    
    q  = q/norm(q);                  % unit quaternion normalization
end 

%% PLOT FIGURES
t       = table(:,1);  
q       = table(:,2:5); 
euler   = rad2deg*table(:,6:8);
w       = rad2deg*table(:,9:11);  
tau     = table(:,12:14);
euler_d = rad2deg*table(:,15:17);
e_tilde = rad2deg*table(:, 18:20);

figure(1);
subplot(3, 1, 1);
    hold on;
    plot(t, euler(:,1), 'b');
    plot(t, euler_d(:,1), 'b--');
    hold off;
    legend('\phi', '\phi_d');
    grid on;
    title('Euler angles');
subplot(3, 1, 2);
    hold on;
    plot(t, euler(:,2), 'r');
    plot(t, euler_d(:,2), 'r--');
    hold off;
    legend('\theta', '\theta_d');
    grid on;
    ylabel('angle [deg]');
subplot(3, 1, 3);
    hold on;
    plot(t, euler(:,3), 'g');
    plot(t, euler_d(:,3), 'g--');
    hold off;
    legend('\psi', '\psi_d');
    grid on;
    xlabel('time [s]'); 


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
brg_plot(t, euler_d);
grid on;
legend('phi', 'theta', 'psi');
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