function result = matlab_abm_run(config, output_directory)
%MATLAB_ABM_RUN Simulate and save one maintained MATLAB ABM run.

arguments
    config (1,1) struct
    output_directory (1,1) string = fullfile("runs", "matlab", config.name)
end

result = matlab_abm_simulate(config);
if ~isfolder(output_directory)
    mkdir(output_directory);
end
save(fullfile(output_directory, 'result.mat'), 'result', '-v7.3');
save(fullfile(output_directory, 'config.resolved.mat'), 'config');
summary_json = jsonencode(result.summary, PrettyPrint=true);
file_id = fopen(fullfile(output_directory, 'summary.json'), 'w');
if file_id < 0
    error('urbancrime:Output', 'Could not create summary.json.');
end
cleanup = onCleanup(@() fclose(file_id));
fprintf(file_id, '%s\n', summary_json);
fprintf('Completed %s; output: %s\n', config.name, output_directory);
end
