%function test_system(option)
%

%if nargin < 1, option = 1; end

option         = 3;
ell            = 1;
delta_t        = 0.01;
omega          = 1/15;
As0            = 1/30;
% parameter setting on page 1257
switch option
  case 1
    eta        = 0.2;
    theta      = 0.56;
    Gamma      = 0.019;
  case 2
    eta        = 0.2;
    theta      = 5.6;
    Gamma      = 0.002;    
  case 3
    eta        = 0.03;
    theta      = 0.56;
    Gamma      = 0.019;    
  case 4
    eta        = 0.03;
    theta      = 5.6;
    Gamma      = 0.002;    
end
% grid size
num_rows       = 128;
num_cols       = 128;
% final time
T_f            = 150;
tspan          = [0, T_f];
% compute the initial conditions
B_bar          = theta * Gamma/omega;
A_bar          = As0 + B_bar;
n_bar          = Gamma * delta_t/(1 - exp(-A_bar * delta_t));
%
params.ell     = ell;
params.delta_t = delta_t;
params.omega   = omega;
params.As0     = As0;
params.eta     = eta;
params.theta   = theta;
params.Gamma   = Gamma;
params.BC_type = 'noFlow';
params.skips   = 100;
% compute the initial condition
B0             = B_bar * ones(num_rows, num_cols);
total_n        = ceil(n_bar * num_rows * num_cols);
% we need to have at least one criminal
if total_n == 0; total_n = 1; end
n_idx          = randi(num_rows * num_cols, total_n);
% set up n0
n0             = zeros(num_rows, num_cols);
for n_count = 1 : total_n
  idx          = n_idx(n_count);
  if n0(idx) == 0
    n0(idx)    = 1;
  else
    n0(idx)    = n0(idx) + 1;
  end
end
% evolve the system
[A, n, E]      = system_evolve(B0, n0, tspan, params);
%
x_min          = 0;
x_max          = num_rows * ell;
y_min          = 0;
y_max          = num_cols * ell;
[X, Y]         = generate_grid(x_min, y_min, num_rows, num_cols, ell);
% plot them in a movie
num_steps      = size(A, 3);
figure('Name', 'Testing', 'Position', [50, 50, 900, 300]);
vidfile        = VideoWriter('CrimeModel.mp4','MPEG-4');
open(vidfile);
for ind = 1:num_steps
    t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;
    [CO, theMap] = remap_colormap(squeeze(A(:, :, ind)), 2);
    surface(X, Y, squeeze(A(:, :, ind)), CO, 'edgeColor', 'none');
    colormap(theMap);
    colorbar
    title('Attractiveness');
    axis tight;
    nexttile;
    [CO, theMap] = remap_colormap(squeeze(E(:, :, ind)), 2);
    surface(X, Y, squeeze(E(:, :, ind)), CO, 'edgeColor', 'none');
    colormap(theMap);
    colorbar
    title('Burglary Events')    
    axis tight;
    nexttile;
    [CO, theMap] = remap_colormap(squeeze(n(:, :, ind)), 2);
    surface(X, Y, squeeze(n(:, :, ind)), CO, 'edgeColor', 'none');
    colormap(theMap);
    colorbar
    title('Number of Agents')
    title(t, sprintf('t = %d secs', ind));
    axis tight;
    drawnow
    theFrame = getframe(gcf); 
    writeVideo(vidfile, theFrame);
end
close(vidfile)
%end