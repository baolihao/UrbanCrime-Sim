function matlab_abm_plot_means(result, output_file)
%MATLAB_ABM_PLOT_MEANS Overlay A, rho, H, and pi spatial averages.

arguments
    result (1,1) struct
    output_file (1,1) string
end

figure_handle = figure('Visible', 'off', 'Color', 'white');
hold on;
plot(result.metrics.time, result.metrics.mean_A, 'LineWidth', 1.2, ...
    'DisplayName', '<A>');
plot(result.metrics.time, result.metrics.mean_rho, 'LineWidth', 1.2, ...
    'DisplayName', '<rho>');
plot(result.metrics.time, result.metrics.mean_H, '--', 'LineWidth', 1.2, ...
    'DisplayName', '<H>');
plot(result.metrics.time, result.metrics.mean_pi, '--', 'LineWidth', 1.2, ...
    'DisplayName', '<pi>');
xlabel('time');
ylabel('spatial average');
grid on;
legend('Location', 'best');
exportgraphics(figure_handle, output_file, 'Resolution', 220);
close(figure_handle);
end
