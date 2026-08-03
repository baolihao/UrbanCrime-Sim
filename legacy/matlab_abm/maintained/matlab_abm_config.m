function config = matlab_abm_config(paper, case_id)
%MATLAB_ABM_CONFIG Paper-scale configuration for the maintained MATLAB ABM.
%   CONFIG = MATLAB_ABM_CONFIG("m3as", CASE_ID) supports M3AS ABM cases 1--3.
%   CONFIG = MATLAB_ABM_CONFIG("police", CASE_ID) supports the extension-paper
%   ABM cases 1, 2, 8, and 9. Values are nondimensional unless they are inside
%   CONFIG.scaling, which only maps integer agents to continuum fields.

arguments
    paper (1,1) string
    case_id (1,1) double {mustBeInteger, mustBePositive}
end

config.schema_version = 2;
config.grid.boundary = "no_flow";
config.time.start = 0.0;
config.time.step = 0.02;
config.time.save_every = 50;
config.scaling.omega = 1/15;
config.scaling.beta = 0.25;
config.model.arrest_probability = 0.0;
config.model.generation_distribution = "bernoulli";
config.initial.seed = 1000;
config.policing.strategy = "none";
config.policing.tau = NaN;

switch lower(paper)
    case "m3as"
        config.name = sprintf('m3as-abm-case-%d', case_id);
        config.grid.nx = 200;
        config.grid.ny = 200;
        config.grid.spacing = 0.08;
        config.scaling.diffusivity = 0.5;
        config.model.ast = 1/30;
        config.model.generation_ratio = 1.0;
        config.initial.dynamic_attractiveness = 1.0;
        config.initial.criminal_density = 0.8;
        config.initial.police_density = 0.0;
        config.initial.information = 0.0;
        config.initial.seed = 1000 + case_id;
        switch case_id
            case 1
                config.time.final = 300.0;
                config.scaling.theta = 0.58;
                config.model.eta = 0.9;
                config.model.generation_ratio = 1.00485;
                config.output.snapshot_times = [0, 3.6, 10, 15, 200, 225, 250, 275, 300];
            case 2
                config.time.final = 200.0;
                config.scaling.theta = 0.2339;
                config.model.eta = 0.3;
                config.model.generation_ratio = 0.9999225;
                config.output.snapshot_times = [0, 5, 10, 20, 200];
            case 3
                config.time.final = 200.0;
                config.scaling.theta = 0.2339;
                config.model.eta = 0.03;
                config.model.generation_ratio = 0.9999225;
                config.output.snapshot_times = [0, 1.6, 10, 20, 200];
            otherwise
                error('urbancrime:UnknownCase', 'M3AS ABM case must be 1, 2, or 3.');
        end

    case {"police", "extension"}
        if ~ismember(case_id, [1, 2, 8, 9])
            error('urbancrime:UnknownCase', ...
                'Extension-paper ABM case must be 1, 2, 8, or 9.');
        end
        config.name = sprintf('police-response-abm-case-%d', case_id);
        config.grid.nx = 100;
        config.grid.ny = 100;
        config.grid.spacing = 0.1;
        config.scaling.diffusivity = 4.0;
        config.model.ast = 0.02;
        config.model.generation_ratio = 1.5;
        config.initial.dynamic_attractiveness = 1.5;
        config.initial.criminal_density = 0.6;
        config.initial.police_density = 0.5;
        config.initial.seed = 2000 + case_id;
        config.policing.strategy = "delayed";
        config.policing.tau = 5.0;
        switch case_id
            case 1
                config.time.final = 150.0;
                config.scaling.theta = 0.58;
                config.model.eta = 0.7;
                config.output.snapshot_times = [0, 150];
            case 2
                config.time.final = 500.0;
                config.scaling.theta = 0.58;
                config.model.eta = 0.3;
                config.output.snapshot_times = [0, 500];
            case 8
                config.time.final = 780.0;
                config.scaling.theta = 0.2339;
                config.model.eta = 0.15;
                config.policing.tau = 50.0;
                config.output.snapshot_times = [0, 765, 768, 771, 774, 777, 780];
            case 9
                config.time.final = 759.0;
                config.scaling.theta = 0.2339;
                config.model.eta = 0.15;
                config.initial.police_density = 0.1;
                config.output.snapshot_times = [0, 744, 747, 750, 753, 756, 759];
        end
        config.initial.information = config.initial.criminal_density .* ...
            (config.model.ast + config.initial.dynamic_attractiveness) .* ...
            exp(-config.initial.police_density);

    otherwise
        error('urbancrime:UnknownPaper', 'Paper must be "m3as" or "police".');
end
end
