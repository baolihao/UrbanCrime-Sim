function result = reproduce_m3as(case_id, output_root)
%REPRODUCE_M3AS Run one of the three M3AS MATLAB ABM cases.

arguments
    case_id (1,1) double {mustBeMember(case_id, [1, 2, 3])} = 1
    output_root (1,1) string = fullfile("runs", "matlab")
end
config = matlab_abm_config("m3as", case_id);
output_directory = fullfile(output_root, config.name);
result = matlab_abm_run(config, output_directory);
switch case_id
    case 1
        matlab_abm_plot_fields(result, [0, 3.6, 10, 15, 200], ...
            ["A", "rho"], fullfile(output_directory, 'figure5.png'));
        matlab_abm_plot_fields(result, [225, 250, 275, 300], ...
            "A", fullfile(output_directory, 'figure6.png'));
    case 2
        matlab_abm_plot_fields(result, [0, 5, 10, 20, 200], ...
            ["A", "rho"], fullfile(output_directory, 'figure9.png'));
    case 3
        matlab_abm_plot_fields(result, [0, 1.6, 10, 20, 200], ...
            ["A", "rho"], fullfile(output_directory, 'figure12.png'));
end
end
