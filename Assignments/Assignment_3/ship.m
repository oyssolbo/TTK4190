function xdot = ship(x,u)
%% Information
% xdot = ship(x,u) returns the time derivative of the state vector: 
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
% Author:           Solbø, Ø. & Strøm, C.
% Study program:    MTTK
% Date:             date

% Check of input and state dimensions
if (length(x)~= 8),error('x-vector must have dimension 8 !');end
if (length(u)~= 2),error('u-vector must have dimension 2 !');end

%% Dimensional states and input
delta_c = u(1); 
n_c     = u(2);

nu    = x(1:3)';
eta   = x(4:6)';
delta = x(7);
n     = x(8); 

%% Ship parameters 
m = 17.0677e6;          % Mass (kg)
Iz = 2.1732e10;         % Yaw moment of inertia (kg m^3)
xg = -3.7;              % CG x-ccordinate (m)
L = 161;                % Length (m)
B = 21.8;               % Beam (m)
T = 8.9;                % Draft (m)
KT = 0.7;               % Propeller coefficient (-)
Dia = 3.3;              % Propeller diameter (m)
rho = 1025;             % Density of water (m/s^3)
T_1 = 20;               % Time constant surge
T_2 = 20;               % Time constant sway
T_6 = 10;               % Time constant yaw

% Rudder limitations
delta_max  = 40 * pi/180;        % Max rudder angle      (rad)
Ddelta_max = 5  * pi/180;        % Max rudder derivative (rad/s)

%% Mass matrix
% Added mass matrix
Xudot = -8.9830e5;
Yvdot = -5.1996e6;
Yrdot =  9.3677e5;
Nvdot =  Yrdot;
Nrdot = -2.4283e10;

M_A = -[Xudot,  0,      0;
        0,      Yvdot,  Yrdot;
        0,      Nvdot,  Nrdot];

% Rigid-body mass matrix
M_RB = [m, 0,    0; 
        0, m,    m*xg;
        0, m*xg, Iz];
% M_inv = inv(M_RB);

% Mass matrix according to 6.103
M = M_A + M_RB;
% M_inv = inv(M);

%% Input matrix
t_thr = 0.05;           % Thrust deduction number
X_delta2 = 0;           % Rudder coefficients (Section 9.5)
Y_delta = 0;      
N_delta = 1;
Bi = [(1-t_thr),  X_delta2;
       0,         Y_delta;
       0,         N_delta];
    
%% Coreolis and centripetal matrix
% Rigid body coreolis matrix
C_RB = m * nu(3) * [0 -1 -xg; 
                    1  0  0;
                    xg 0  0];

% Hydrodynamic coreolis matrix for ALL states
C_A_11 = zeros(3);
C_A_12 = Smtrx([Xudot*nu(1); Yvdot*nu(2)+Yrdot*nu(3); 0]);
C_A_21 = C_A_12;
C_A_22 = Smtrx([0; 0; Nvdot*nu(2) + Nrdot*nu(3)]);

% Hydrodynamic corelois matrix for DOF 1,2,6
C_A = [C_A_11(1:2,1:2), C_A_12(1:2,3);
       C_A_21(3,1:2),   C_A_22(3,3)];

% Coreolis and centripetal matrix according to 6.104
C = C_A + C_RB;

%% Damping
% Linear damping
Xu = -(m - Xudot)/T_1;
Yv = -(m - Yvdot)/T_2;
Nr = -(Iz - Nrdot)/T_6;

D = [Xu, 0,  0;
     0,  Yv, 0;
     0,  0,  Nr];

% Nonlinear drag in surge
S = 2*T*L + 2*B*T + B*L; % Wetted surface of the hull (sides + keel)
k = 0.1;
epsilon = 0.001;
C_R = 0;

L_pp = L; % Characteristic linear dimension (length of the ship) 
% abs(nu(1)) gets fucked when surge speed explodes
R_n = L_pp*1e6*abs(nu(1));
C_f = 0.075/((log10(R_n) - 2)^2 + epsilon) + C_R;

D_X = -1/2*rho*S*(1+k)*C_f*abs(nu(1))*nu(1);

% Cross-flow drag
num_strips = 10;
dx = L/num_strips; 
Cd_2D = Hoerner(B,T);

Ycf = 0; % Initial cross-flow drag in Y
Ncf = 0; % Initial cross-flow drag in N

v_r = nu(2);
for xL = -L/2:dx:L/2,
    Ucf = abs(v_r + xL*nu(3)) * (v_r + xL*nu(3));
    Ycf = Ycf - 0.5*rho*T*Cd_2D*Ucf*dx;    % Sway force
    Ncf = Ncf - 0.5*rho*T*Cd_2D*xL*Ucf*dx; % Yaw moment
end

D_n = [D_X, 0,   0; 
       0,   Ycf, 0;
       0,   0,   Ncf];

% Total drag
D = D + D_n;

%% Rotation matrix
R = Rzyx(0,0,eta(3));

%% Thrust 
thr = rho * Dia^4 * KT * abs(n) * n;    % Thrust command (N)

%% Ship dynamics
u = [ thr delta ]';
tau = Bi * u;
nu_dot = M \ (tau - C*nu - D*nu); % M_inv * (tau - C_RB * nu); 
eta_dot = R * nu;    

% Rudder saturation and dynamics (Sections 9.5.2)
if abs(delta_c) >= delta_max
    delta_c = sign(delta_c)*delta_max;
end

delta_dot = delta_c - delta;
if abs(delta_dot) >= Ddelta_max
    delta_dot = sign(delta_dot)*Ddelta_max;
end    

%% Propeller dynamics
n_dot = (1/10) * (n_c - n);

%% Return value
xdot = [nu_dot' eta_dot' delta_dot n_dot]';
end