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

% Dimensional states and input
delta_c = u(1); 
n_c     = u(2);

nu    = x(1:3)';
eta   = x(4:6)';
delta = x(7);
n     = x(8); 

% ship parameters 
m = 17.0677e6;          % mass (kg)
Iz = 2.1732e10;         % yaw moment of inertia (kg m^3)
xg = -3.7;              % CG x-ccordinate (m)
L = 161;                % length (m)
B = 21.8;               % beam (m)
T = 8.9;                % draft (m)
KT = 0.7;               % propeller coefficient (-)
Dia = 3.3;              % propeller diameter (m)
rho = 1025;             % density of water (m/s^3)

% rudder limitations
delta_max  = 40 * pi/180;        % max rudder angle      (rad)
Ddelta_max = 5  * pi/180;        % max rudder derivative (rad/s)

% added mass matrix
Xudot = -8.9830e5;
Yvdot = -5.1996e6;
Yrdot =  9.3677e5;
Nvdot =  Yrdot;
Nrdot = -2.4283e10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add added mass here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_A = [Xudot       0     0;
           0   Yvdot Yrdot;
           0   Nvdot Nrdot;];

% rigid-body mass matrix
MRB = [ m 0    0 
        0 m    m*xg
        0 m*xg Iz ];

M = MRB + M_A;
Minv = inv(M);

% input matrix
t_thr = 0.05;           % thrust deduction number
X_delta2 = 0;           % rudder coefficients (Section 9.5)
Y_delta = 0;      
N_delta = 1;
Bi = [ (1-t_thr)  X_delta2
        0        Y_delta
        0        N_delta  ];
    
% state-dependent time-varying matrices
CRB = m * nu(3) * [ 0 -1 -xg 
                    1  0  0 
                    xg 0  0  ];
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add Coriolis due to added mass here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a1 = Xudot*x(1);
a2 = Yvdot*x(2);

C_A = [  0  0  a2;
         0  0 -a1;
       -a2 a1   0;];
   
C = CRB + C_A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add linear damping here (D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T1 = 20;
T2 = 20;
T6 = 10;
Xu = -(m - Xudot)/T1;
Yv = -(m - Yvdot)/T2;
Nr = -(Iz - Nrdot)/T6;

D = [Xu   0  0;
      0  Yv  0;
      0   0 Nr;];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add nonlinear damping here (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Surge damping
u_c = 0;
u_r = x(1) - u_c;
nu_c = 1e-6;
epsilon = 0.001;
k = 0.1;
S = L*B + 2*B*T + 2*L*T;

R_n = @(u) L/nu_c *abs(u);
C_f = 0.075/((log10(R_n(u_r)) - 2)^2 + epsilon);

Xn = -0.5*rho*S*(1+k)*C_f*abs(u_r)*u_r;

% Cross-flow drag
v_c = 0;
v_r = x(2) - v_c;
r = x(3);
Cd_2D = Hoerner(B, T);

% integrals from problem text.
dx = L/10;
Ycf = 0;
Ncf = 0;
for xL = -L/2:dx:L/2
    Ucf = abs(v_r + xL * r) * (v_r + xL * r);
    Ycf = Ycf - 0.5 * rho * T * Cd_2D * Ucf * dx; % sway force
    Ncf = Ncf - 0.5 * rho * T * Cd_2D * xL * Ucf * dx; % yaw moment
end

% Final nonlinear damping
Dn = -diag([Xn -Ycf -Ncf]);

R = Rzyx(0,0,eta(3));

% thrust 
thr = rho * Dia^4 * KT * abs(n) * n;    % thrust command (N)

% ship dynamics
u = [ thr delta ]';
tau = Bi * u;
nu_dot = Minv * (tau - C * nu - D*nu - Dn*nu); 
eta_dot = R * nu;    

% Rudder saturation and dynamics (Sections 9.5.2)
if abs(delta_c) >= delta_max
    delta_c = sign(delta_c)*delta_max;
end

delta_dot = delta_c - delta;
if abs(delta_dot) >= Ddelta_max
    delta_dot = sign(delta_dot)*Ddelta_max;
end    

% propeller dynamics
n_dot = (1/10) * (n_c - n);

xdot = [nu_dot' eta_dot' delta_dot n_dot]';

end