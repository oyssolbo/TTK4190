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

%% USER INPUTS
h = 0.1;                     % sample time (s)
N  = 5000;                    % number of samples. Should be adjusted

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
A_e2 = 15;
A_e3 = 10;
T_e2 = 0.1;
T_e3 = 0.05;
t = 1:N + 1;
eps1_d = 0 * t * deg2rad;
eps2_d = A_e2 * cos(T_e2 * t) * deg2rad;
eps3_d = A_e3 * sin(T_e3 * t) * deg2rad;

eta_d = sqrt(1 - eps1_d.^2 - eps2_d.^2 - eps3_d.^2);

q_d = [eta_d; eps1_d; eps2_d; eps3_d]';
q_d_conj = quatconj(q_d);

w_d = [0, 0, 0]';

%% FOR-END LOOP
for i = 1:N+1
   t = (i-1)*h;                  % time
   
   [phi,theta,psi] = q2euler(q); % transform q to Euler angles
   [J,J1,J2] = quatern(q);       % kinematic transformation matrices
   
   % Control law for problem 1.2
   %e = q(2:4);
   %tau = - K_d * w - k_p * e;
   
   % Control law for problem 1.5 and 1.6 (quaternion)
   q_tilde = quatmultiply(q_d_conj(i,:), q');
   e_tilde = q_tilde(2:4);
   %tau = -K_d * w - k_p * e_tilde';
   
   % Control law for problem 1.6 (angular velocity)
   T_euler_inv = [1, 0, -sin(theta);
                  0, cos(phi), cos(theta) * sin(phi);
                  0, -sin(phi), cos(theta) * cos(phi)];

   % Since euler_ang is defined explicitly and in simple closed form, we
   % can manually time differentiate each term to get euler_ang_d without
   % much effort.
   dTheta_d = [0,  -T_e2*A_e2*sin(T_e2*t), T_e3*A_e3*cos(T_e3*t)];
   
   w_d = T_euler_inv * dTheta_d';
   tau = -K_d * (w - w_d) - k_p * e_tilde';
   
   tracking_error(i,:) = e_tilde;
   ref_signal(i,:) = dTheta_d;
   
   q_dot = J2*w;                        % quaternion kinematics
   w_dot = I_inv*(Smtrx(I*w)*w + tau);  % rigid-body kinetics
   
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


figure (1); clf;
hold on;
plot(t, phi, 'b');
plot(t, theta, 'r');
plot(t, psi, 'g');
hold off;
grid on;
legend('\phi', '\theta', '\psi');
title('Euler angles');
xlabel('time [s]'); 
ylabel('angle [deg]');

figure (2); clf;
hold on;
plot(t, w(:,1), 'b');
plot(t, w(:,2), 'r');
plot(t, w(:,3), 'g');
hold off;
grid on;
legend('p', 'q', 'r');
title('Angular velocities');
xlabel('time [s]'); 
ylabel('angular rate [deg/s]');

figure (3); clf;
hold on;
plot(t, tau(:,1), 'b');
plot(t, tau(:,2), 'r');
plot(t, tau(:,3), 'g');
hold off;
grid on;
legend('x', 'y', 'z');
title('Control input');
xlabel('time [s]'); 
ylabel('input [Nm]');

figure (4); clf;
hold on;
plot(t, e_tilde(:,1), 'b');
plot(t, e_tilde(:,2), 'r');
plot(t, e_tilde(:,3), 'g');
hold off;
grid on;
legend('e_1', 'e_2', 'e_3');
title('Tracking error - epsilon');
xlabel('time [s]'); 
ylabel('angle [deg]');

figure (5); clf;
hold on;
plot(t, ref_signal(:,1), 'b');
plot(t, ref_signal(:,2), 'r');
plot(t, ref_signal(:,3), 'g');
hold off;
grid on;
legend('e_1', 'e_2', 'e_3');
title('Reference signal');
xlabel('time [s]'); 
ylabel('angle [deg]');