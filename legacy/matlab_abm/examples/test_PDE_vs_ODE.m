% default
option            = 2;
% the space domain is fixed at [0, L]x[0, L], L_tilde = 16 (this is the total length)
L_tilde           = 16;
x_min_tilde       = 0;
x_max_tilde       = L_tilde;
y_min_tilde       = 0;
y_max_tilde       = L_tilde;
% the time domain is fixed at [0, T], T = 200
T_tilde           = 200;
% choose the spatial mesh step size, and time step size
num_rows          = 200;
num_cols          = num_rows;
ell_tilde         = x_max_tilde/num_rows;
delta_t_tilde     = 1/50;
% these three cases are from the note on overleaf
% again (**quan**)_tilde means dimensionless variables
switch option
  case 1
% no hotspots    
    eta           = 0.9;
    Ast_tilde     = 1/30;
    % B0_tilde = Gamma * theta/omega^2
    B0_tilde      = 1.0;
    % A0_tilde = B0_tilde + Ast_tilde
    A0_tilde      = (B0_tilde + Ast_tilde) * ones(num_rows, num_cols);
    ind           = A0_tilde < 0;
    A0_tilde(ind) = 0;    
    % rho0 = 0.8
    rho0_tilde     = 0.8 * ones(num_rows, num_cols);
    id             = rho0_tilde < 0;
    rho0_tilde(id) = 0;
  case 2
% small hotspots    
    eta           = 0.3;
    Ast_tilde     = 1/30;
    % B0_tilde = Gamma * theta/omega^2
    B0_tilde      = 1.0;
    % A0_tilde = B0_tilde + Ast_tilde + random noise
    A0_tilde      = B0_tilde + Ast_tilde + normrnd(0, 0.05, num_rows, num_cols);
    ind           = A0_tilde < 0;
    A0_tilde(ind) = 0;    
    % rho0 = 0.6 + random noise
    rho0_tilde     = 0.8 + normrnd(0, 0.01, num_rows, num_cols);
    id             = rho0_tilde < 0;
    rho0_tilde(id) = 0;
  case 3
% big hotspots    
    eta           = 0.03;
    Ast_tilde     = 1/30;
    % B0_tilde = Gamma * theta/omega^2
    B0_tilde      = 1.0;
    % A0_tilde = B0_tilde + Ast_tilde + random noise
    A0_tilde      = B0_tilde + Ast_tilde + normrnd(0, 0.05^2, num_rows, num_cols);
    ind           = A0_tilde < 0;
    A0_tilde(ind) = 0;    
    % rho0 = 0.6 + random noise
    rho0_tilde     = 0.8 + normrnd(0, 0.01, num_rows, num_cols);
    id             = rho0_tilde < 0;
    rho0_tilde(id) = 0;
  case 4
% this is the case 2 but with Ast and B0 being maps    
    eta           = 0.3;
    [X_tilde, ...
      Y_tilde]    = generate_grid(x_min_tilde, y_min_tilde, num_rows, num_cols, ell_tilde);
% sparsity level in percentage
% here we take sparsity = 1 - delta in the note
    sparsity      = 0;
    Ast_tilde     = 1/30 * mask_from_sparsity(sparsity, ones(num_rows, num_cols)) ...
                    + 1/30 * get_highway(X_tilde, Y_tilde, L_tilde);
    % B0_tilde = Gamma * theta/omega^2
    sparsity      = 0.9;
    B0_tilde      = mask_from_sparsity(sparsity, ones(num_rows, num_cols)) ...
                    + get_highway(X_tilde, Y_tilde, L_tilde); 
    % A0_tilde = B0_tilde + Ast_tilde + random noise
    A0_tilde      = B0_tilde + Ast_tilde + normrnd(0, 0.1, num_rows, num_cols);
    ind           = A0_tilde < 0;
    A0_tilde(ind) = 0;    
    % rho0 = 0.6 + random noise
    rho0_tilde     = 0.8 + normrnd(0, 0.01, num_rows, num_cols);
    id             = rho0_tilde < 0;
    rho0_tilde(id) = 0;
  otherwise
end
%
[A_tilde, n_tilde, E, rho_tilde, ts_tilde, X_tilde, Y_tilde]  = compare_system(Ast_tilde, ...
  A0_tilde, rho0_tilde, eta, B0_tilde, delta_t_tilde, ell_tilde, T_tilde, option);
% plot them in a movie
num_steps          = size(A_tilde, 3);
figure('Name', 'Dimension Output', 'Position', [50, 50, 900, 300]);
vidfile            = VideoWriter('CrimeModel.mp4','MPEG-4');
open(vidfile);
% when we plot, we change them back to dimensionless
for ind = 1 : num_steps
    TL = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;
    [CO, theMap] = remap_colormap(squeeze(A_tilde(:, :, ind)), 2);
    surface(X_tilde, Y_tilde, squeeze(A_tilde(:, :, ind)), CO, 'edgeColor', 'none');
    colormap(theMap)
    colorbar
    title('Attractiveness')
    axis tight;
    nexttile;
    [CO, theMap] = remap_colormap(squeeze(E(:, :, ind)));
    surface(X_tilde, Y_tilde, squeeze(E(:, :, ind)), CO, 'edgeColor', 'none');
    colormap(theMap)
    colorbar
    title('Burglary Events')
    axis tight;
    nexttile;
    [CO, theMap] = remap_colormap(squeeze(n_tilde(:, :, ind)));
    surface(X_tilde, Y_tilde, squeeze(n_tilde(:, :, ind)), CO, 'edgeColor', 'none');
    colormap(theMap)
    colorbar
    title('Number of Agents')
    axis tight;
    nexttile;
    [CO, theMap] = remap_colormap(squeeze(rho_tilde(:, :, ind)));
    surface(X_tilde, Y_tilde, squeeze(rho_tilde(:, :, ind)), CO, 'edgeColor', 'none');
    colormap(theMap)
    colorbar
    title('Agent Density')
    axis tight;    
    title(TL, sprintf('t = %d secs', ts_tilde(ind)));
    drawnow
    theFrame   = getframe(gcf); 
    writeVideo(vidfile, theFrame);
end
close(vidfile)
% end if this is a function instead of a script