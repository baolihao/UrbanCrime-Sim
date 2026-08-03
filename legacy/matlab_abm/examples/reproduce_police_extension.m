function result = reproduce_police_extension(case_id, output_root)
%REPRODUCE_POLICE_EXTENSION Run an extension-paper MATLAB ABM case.

arguments
    case_id (1,1) double {mustBeMember(case_id, [1, 2, 8, 9])} = 1
    output_root (1,1) string = fullfile("runs", "matlab")
end
config = matlab_abm_config("police", case_id);
output_directory = fullfile(output_root, config.name);
result = matlab_abm_run(config, output_directory);
switch case_id
    case 1
        matlab_abm_plot_fields(result, 150, ["A", "H", "rho", "pi"], ...
            fullfile(output_directory, 'figure5.png'));
    case 2
        matlab_abm_plot_means(result, fullfile(output_directory, 'figure6.png'));
    case 8
        matlab_abm_plot_fields(result, [765, 768, 771, 774, 777, 780], ...
            "expected_S", fullfile(output_directory, 'figure18.png'));
    case 9
        matlab_abm_plot_fields(result, [744, 747, 750, 753, 756, 759], ...
            "expected_S", fullfile(output_directory, 'figure22.png'));
end
end
