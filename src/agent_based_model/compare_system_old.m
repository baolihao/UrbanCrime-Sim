function compare_system(As0, A0, rho0, eta, B_bar, delta_t, ell, T_f, option)
%

% here As0 is the A0 in the paper, A0 is the initial of A, rho0 is the initial of rho
% b_bar = Gamma * theta/omega^2

if nargin < 9, option = 1; end
% trying to keep the same CFD rate
D              = ell^2/delta_t;
omega          = 1/15;
z              = 4;
num_rows       = size(A0, 1);
num_cols       = size(A0, 2);
% parameter setting on page 1257
switch option
  case 1
    theta      = 0.56;
  case 2
    theta      = 5.6;    
end
% switch it back to dimension variables
delta_t        = delta_t/omega;
ell            = sqrt(D/(z * omega)) * ell;
epsilon        = theta * delta_t;
Gamma          = omega^2/theta * B_bar;
% final time
T_f            = T_f/omega;
tspan          = [0, T_f];
x_min          = 0;
%x_max          = num_rows * ell;
y_min          = 0;
%y_max          = num_cols * ell;
% compute the initial conditions
A0             = A0 * omega;
B0             = A0 - As0 * omega;
rho0           = rho0 * omega/(epsilon * D);
n0             = ceil(rho0 * ell^2);
% pack the parameters
params.ell     = ell;
params.delta_t = delta_t;
params.omega   = omega;
params.As0     = As0;
params.eta     = eta;
params.theta   = theta;
params.Gamma   = Gamma;
params.BC_type = 'noFlow';
params.skips   = 150;
% evolve the system
[A, n, E]      = system_evolve(B0, n0, tspan, params);
[X, Y]         = generate_grid(x_min, y_min, num_rows, num_cols, ell);
% plot them in a movie
num_steps      = size(A, 3);
figure('Name', 'Dimension Output', 'Position', [50, 50, 900, 300]);
vidfile        = VideoWriter('CrimeModel.mp4','MPEG-4');
open(vidfile);
for ind = 1 : num_steps
    t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;
    [CO, theMap] = remap_colormap(squeeze(A(:, :, ind)), 2);
    surface(X, Y, squeeze(A(:, :, ind)), CO, 'edgeColor', 'none');
    colormap(theMap)
    colorbar
    title('Attractiveness')
    axis tight;
    nexttile;
    [CO, theMap] = remap_colormap(squeeze(E(:, :, ind)));
    surface(X, Y, squeeze(E(:, :, ind)), CO, 'edgeColor', 'none');
    colormap(theMap)
    colorbar
    title('Burglary Events')
    axis tight;
    nexttile;
    [CO, theMap] = remap_colormap(squeeze(n(:, :, ind)));
    surface(X, Y, squeeze(n(:, :, ind)), CO, 'edgeColor', 'none');
    colormap(theMap)
    colorbar
    title('Number of Agents')
    axis tight;
    title(t, sprintf('t = %d secs', ind));
    drawnow
    theFrame   = getframe(gcf); 
    writeVideo(vidfile, theFrame);
end
close(vidfile)
end