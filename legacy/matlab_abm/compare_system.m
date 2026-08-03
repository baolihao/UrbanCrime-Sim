function [A_tilde, n_tilde, E, rho_tilde, ts_tilde, X_tilde, Y_tilde] = compare_system(Ast_tilde, ...
          A0_tilde, rho0_tilde, eta, B0_tilde, delta_t_tilde, ell_tilde, T_tilde, option)
% function compare_system(Ast_tilde, A0_tilde, rho0_tilde, eta, B_bar_tilde, delta_t_tilde, 
% ell_tilde, T_tilde, option)
% (**quant**)_tilde ==> dimensionless variables from PDE
% Here 
% Ast_tilde: the static (intrinsic) part of A_tilde
% A0_tilde: initial (dimensionless) attractivenes
% rho0_tilde: initial (dimensionless) density
% eta: measure neighborhood effects, eta \in [0, 1] 
% B0_tilde = B0 * omega, where B0 = Gamma * theta/omega
% delta_t_tilde: (dimensionless) time step size
% ell_tilde: (dimensionless) spatial mesh step size
% T_tilde: (dimensionless) final time
% option: choosing theta
% we need to choose omega, theta, delta_t, and D in order to run the ODE code

if nargin < 9, option = 1; end
% choosing D, omega, and theta are very important, they should be set first
% D = ell^2/delta_t (reciprical of the CFL number)
D                  = 0.5;
% dynamic attractiveness decay rate
omega              = 1/15;
% parameter setting on page 1257
switch option
% theta: increase in attractiveness due to one burglary event  
  case 1
    theta          = 0.56;
  case 2
    theta          = 5.6;    
end
% number of neighgors in a 2D grid
z                  = 4;
% time step size, better to be just delta_t = delta_t_tilde/omega
delta_t            = delta_t_tilde/omega;
% resolution
num_rows           = size(A0_tilde, 1);
num_cols           = size(A0_tilde, 2);
% characteristic length
ell_c              = sqrt(D/omega);
% controlling time step size
epsilon            = theta * delta_t;
% scaling factor to change spatial mesh to dimensional less
x_scale            = sqrt(z)/ell_c;
rho_scale          = epsilon * ell_c^2;
% spatial mesh step size in dimension case
ell                = ell_tilde/x_scale;
% rate of burglar generation at each site
Gamma              = omega^2/theta * B0_tilde;
% final time
T                  = T_tilde/omega;
tspan              = [0, T];
x_min              = 0;
y_min              = 0;
% compute the initial conditions
A0                 = A0_tilde * omega;
Ast                = Ast_tilde *  omega;
B0                 = A0 - Ast;
rho0               = rho0_tilde/rho_scale;
% using ceil to make sure n0 has only integers
n0                 = ceil(rho0 * ell^2);
% packakge the parameters
params.ell         = ell;
params.delta_t     = delta_t;
params.omega       = omega;
params.Ast         = Ast;
params.eta         = eta;
params.theta       = theta;
params.Gamma       = Gamma;
params.BC_type     = 'noFlow';
params.skips       = 100;
% this indicator for whether Ast and Gamma are maps
if isscalar(Ast) && isscalar(Gamma)
  params.use_map   = false;
else
  params.use_map   = true;
end
% evolve the system
[A, n, E, ts]      = system_evolve(B0, n0, tspan, params);
[X, Y]             = generate_grid(x_min, y_min, num_rows, num_cols, ell);
% compute rho
rho                = n/ell^2;
% convert them back to dimensionless
% see Eqn (3.8) on page 1259 of the paper
X_tilde            = X * x_scale;
Y_tilde            = Y * x_scale;
A_tilde            = A/omega;
rho_tilde          = rho * rho_scale;
ts_tilde           = ts * omega;
% this is done by definition of n and rho
n_tilde            = rho_tilde * ell_tilde;
end