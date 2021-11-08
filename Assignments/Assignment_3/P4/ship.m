function [xdot,u] = ship(x,u,nu_c,tau_ext, U_ref)
%% Information
% [xdot,u] = ship(x,u,nu_c,tau_ext) returns the time derivative of the state vector: 
% x = [ u v r x y psi delta n ]' for a ship with L = 161 m where:
%
% u     = surge velocity, must be positive  (m/s)    
% v     = sway velocity                     (m/s)
% r     = yaw velocity                      (rad/s)
% x     = position in x-direction           (m)
% y     = position in y-direction           (m)
% psi   = yaw angle                         (rad)
% delta = actual rudder angle               (rad)
% n     = actual shaft velocity             (rpm)
% 
% The input vector is :
%
% u       = [ delta_c  n_c ]'  where
%
% delta_c = commanded rudder angle          (rad)
% n_c     = commanded shaft velocity        (rpm)
%
% The current vector is : 
%
% nu_c    = [ u_c v_c 0 ]'  where
%
% uc = current velocity in surge            (m/s)
% vc = current velocity in sway             (m/s)
%
% The external environmental force vector is denoted by tau_ext
%
% Author:           Solbø, Ø. & Strøm, C.
% Study program:    MTTK
% Date:             October 2021

%% Assertion
% Check of input and state dimensions
if (length(x)~= 8),error('x-vector must have dimension 8 !');end
if (length(u)~= 2),error('u-vector must have dimension 2 !');end

%% Dimensional states and input
delta_c = u(1); 
n_c     = u(2);

nu    = x(1:3);
eta   = x(4:6);
delta = x(7);
n     = x(8); 

nu_r = nu - nu_c;
uc = nu_c(1);
vc = nu_c(2);

%% Ship parameters 
m = 17.0677e6;          % Mass (kg)
Iz = 2.1732e10;         % Yaw moment of inertia (kg m^2)
xg = -3.7;              % CG x-ccordinate (m)
L = 161;                % Length (m)
B = 21.8;               % Beam (m)
T = 8.9;                % Draft (m)
% KT = 0.7;               % Propeller coefficient (-)
Dia = 3.3;              % Propeller diameter (m)
rho = 1025;             % Density of water (m/s^3)

% Rudder limitations
delta_max  = 40 * pi/180;        % Max rudder angle      (rad)
Ddelta_max = 5  * pi/180;        % Max rudder derivative (rad/s)

% Propeller - cofficients
PD = 1.5;   % Pitch / diameter ratio
AEAO = 0.65;% Blade area ratio
z = 4;      % Num propeller-blades
[K_T, K_Q] = wageningen(0, PD, AEAO, z); 

%% Mass matrix
% Added mass matrix
Xudot = -8.9830e5;
Yvdot = -5.1996e6;
Yrdot =  9.3677e5;
Nvdot =  Yrdot;
Nrdot = -2.4283e10;
M_A = -[Xudot,  0,     0; 
        0,     Yvdot, Yrdot; 
        0,     Nvdot, Nrdot];

% Rigid-body mass matrix
M_RB = [m,  0,    0; 
        0,  m,   m*xg;
        0, m*xg,  Iz];
   
% Total mass matrix
M = M_RB + M_A;

%% Rudder coefficients
b = 2;
AR = 8;
CB = 0.8;

lambda = b^2 / AR;
tR = 0.45 - 0.28*CB;
CN = 6.13*lambda / (lambda + 2.25);
aH = 0.75;
xH = -0.4 * L;
xR = -0.5 * L;

%% Input matrix
t_thr = 0.05;                                        % Thrust deduction number
X_delta2 = 0.5 * (1 - tR) * rho * AR * CN;           % Rudder coefficients (Section 9.5)
Y_delta = 0.25 * (1 + aH) * rho * AR * CN; 
N_delta = 0.25 * (xR + aH*xH) * rho * AR * CN;   

B_i = @(u_r,delta) [(1-t_thr), -u_r^2*X_delta2*delta;
                      0,       -u_r^2*Y_delta;
                      0,       -u_r^2*N_delta];
   
%% Centripetal Coreolis-matrices
% Rigid-body Coreolis-matrix
C_RB = m * nu_r(3) * [0, -1, -xg; 
                      1,  0,  0; 
                      xg, 0,  0];

% Added mass Coreolis-matrix
C_A = [             0,                   0,         Yvdot*nu_r(2) + Yrdot*nu_r(3); 
                    0,                   0,                 -Xudot*nu_r(1); 
      -Yvdot*nu_r(2)-Yrdot*nu_r(3), Xudot*nu_r(1),                0];

%% Damping
% Linear damping
T1 = 20; T2 = 20; T6 = 10;
Xu = -(m-Xudot)/T1;
Yv = -(m-Yvdot)/T2;
Nr = -(Iz-Nrdot)/T6;
D = -diag([Xu Yv Nr]);

% Nonlinear surge damping
eps = 0.001;
CR = 0;
k = 0.1;
S = B*L + 2*T*(B+L);
v = 1e-6;
Rn = L / v * abs(nu_r(1));
Cf = 0.075 / (log10(Rn) - 2 + eps)^2 + CR;
Xns = -0.5*rho*S*(1+k)*Cf*abs(nu_r(1))*nu_r(1);

% Nonlinear cross-flow drag
Cd_2d = Hoerner(B,T);
dx = L/10;
Ycf = 0; Ncf = 0;
for xL = -L/2:dx:L/2
    Ucf = abs(nu_r(2) + xL * nu_r(3)) * (nu_r(2) + xL * nu_r(3));
    Ycf = Ycf - 0.5 * rho * T * Cd_2d * Ucf * dx;
    Ncf = Ncf - 0.5 * rho * T * Cd_2d * xL * Ucf * dx;
end
d = -[Xns Ycf Ncf]';

%% Linearized model using sway-yaw subsystem
% Linear matrices
C_RB_star = zeros(3);
C_RB_star(2,3) = m*U_ref;
C_RB_star(3,3) = m*U_ref*xg;

C_A_star = zeros(3);
C_A_star(3,2) = -Yvdot*U_ref + Xudot*U_ref;
C_A_star(2,3) = -Xudot*U_ref;
C_A_star(3,3) = -Yrdot*U_ref;

N = C_RB_star + C_A_star + D;
b = 2*U_ref*[-Y_delta, -N_delta]';

% Reducing down to sway-yaw subsystem
N_lin = N(2:3,2:3);
M_lin = M(2:3,2:3);

% Transfer function
[den, num] = ss2tf(-M_lin\N_lin, M_lin\b, [0, 1], 0);

%% Propultion system
I_m = 1e5;
K_m = 0.6;
T_m = 10;
tau = 0;
rho = 1025;
epsilon = 1e-14;

w = 0.25;   % Guessed value (0.2 - 0.4) normally
t = 0.1;    % Guessed value (0.05 - 0.2) normally

J_a = @(n_k, u_a) u_a / (n_k*Dia + epsilon);  % Adding epsilong to prevent division by zero
T = @(n_k, u_a) rho*Dia^4*K_T*J_a(n_k, u_a)*abs(n_k)*n_k;
Q = @(n_k, u_a) rho*Dia^5*K_Q*J_a(n_k, u_a)*abs(n_k)*n_k;

u_a = (1 - w) * x(1);

% Q = K_Q * rho * D^5 * abs(n) * n;
% T = K_T * rho * D^4 * abs(n) * n;
T_d = U_ref*Xu / (t - 1);

Q_d = @(n_d, u_a) Q(n_d, u_a); 

n_d = sqrt(abs(T_d / (K_T*rho*Dia^4)))*sign(T_d);
Q_m = 1/(T_m + 1) * Q_d(n_d, u_a);

n_dot = 1/I_m * (Q_m - Q(n, u_a)); % Converting n from rpm to rps?

%thr = T(n, u_a); %rho * Dia^4 * K_T * abs(n) * n;    % Thrust command (N)
thr = T_d;
%% Ship dynamics
R = Rzyx(0,0,eta(3));
u = [ thr delta ]';
tau = B_i(nu_r(1),delta) * u;
nu_dot = [nu(3)*vc -nu(3)*uc 0]' + M \ (tau_ext + tau - (C_RB + C_A + D) * nu_r - d);
eta_dot = R * nu;    

% Rudder saturation and dynamics (Sections 9.5.2)
if abs(delta_c) >= delta_max
    delta_c = sign(delta_c)*delta_max;
end

delta_dot = delta_c - delta;
if abs(delta_dot) >= Ddelta_max
    delta_dot = sign(delta_dot)*Ddelta_max;
end    

% Propeller dynamics
% n_dot = (1/10) * (n_c - n);

xdot = [nu_dot' eta_dot' delta_dot n_dot]';
u = [delta_c n_c]';
end