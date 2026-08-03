function result = matlab_abm_simulate(config)
%MATLAB_ABM_SIMULATE Discrete stochastic burglary ABM in paper variables.
% Criminals and police are integer site populations. Saved A, rho, pi, H, and
% S fields use the nondimensional variables printed in the two papers.

validate_config(config);
rng(config.initial.seed, 'twister');

ny = config.grid.ny;
nx = config.grid.nx;
h2 = config.grid.spacing^2;
dt = config.time.step;
omega = config.scaling.omega;
theta = config.scaling.theta;
beta = config.scaling.beta;
diffusivity = config.scaling.diffusivity;

ast = expand_field(config.model.ast, ny, nx, 'model.ast');
generation_ratio = expand_field( ...
    config.model.generation_ratio, ny, nx, 'model.generation_ratio');
B = expand_field(config.initial.dynamic_attractiveness, ny, nx, ...
    'initial.dynamic_attractiveness');
rho0 = expand_field(config.initial.criminal_density, ny, nx, ...
    'initial.criminal_density');
pi0 = expand_field(config.initial.police_density, ny, nx, ...
    'initial.police_density');

rho_unit = 4 * theta * dt / (omega * h2);
pi_unit = 4 * beta * omega / (diffusivity * h2);
event_unit = 4 * theta / (omega * h2);
event_increment = theta / omega;
generation_scale = omega / theta;

criminals = matlab_abm_integer_population(rho0 / rho_unit);
if config.policing.strategy == "none"
    if any(pi0 ~= 0, 'all')
        error('urbancrime:InvalidInitialPolice', ...
            'No-police configurations require zero initial police density.');
    end
    police = zeros(ny, nx);
    H = zeros(ny, nx);
else
    police = matlab_abm_integer_population(pi0 / pi_unit);
    H = expand_field(config.initial.information, ny, nx, 'initial.information');
end

num_steps = round((config.time.final - config.time.start) / dt);
snapshot_steps = round((config.output.snapshot_times - config.time.start) / dt);
snapshot_steps = unique([0, snapshot_steps, num_steps]);
snapshot_slot = zeros(1, num_steps + 1);
snapshot_slot(snapshot_steps + 1) = 1:numel(snapshot_steps);
metric_steps = unique([0:config.time.save_every:num_steps, num_steps]);
metric_slot = zeros(1, num_steps + 1);
metric_slot(metric_steps + 1) = 1:numel(metric_steps);

ns = numel(snapshot_steps);
nm = numel(metric_steps);
result.time = config.time.start + snapshot_steps * dt;
result.A = zeros(ny, nx, ns);
result.rho = zeros(ny, nx, ns);
result.events = zeros(ny, nx, ns);
result.pi = zeros(ny, nx, ns);
result.H = zeros(ny, nx, ns);
result.expected_S = zeros(ny, nx, ns);
result.realized_S = zeros(ny, nx, ns);
result.criminal_counts = zeros(ny, nx, ns);
result.police_counts = zeros(ny, nx, ns);
result.metrics.time = config.time.start + metric_steps * dt;
result.metrics.mean_A = zeros(1, nm);
result.metrics.mean_rho = zeros(1, nm);
result.metrics.mean_pi = zeros(1, nm);
result.metrics.mean_H = zeros(1, nm);
result.metrics.mean_expected_S = zeros(1, nm);
result.metrics.mean_realized_S = zeros(1, nm);
result.metrics.cumulative_events = zeros(1, nm);

last_events = zeros(ny, nx);
last_realized = zeros(ny, nx);
total_events = 0;
total_generated = 0;
initial_criminals = sum(criminals, 'all');
initial_police = sum(police, 'all');
save_snapshot(1);
save_metrics(1);

for step = 1:num_steps
    A = ast + B;
    pi = pi_unit * police;
    probability = -expm1(-A .* exp(-pi) * dt);
    probability = min(max(probability, 0), 1);
    events = binornd(criminals, probability);
    criminal_arrivals = move_population( ...
        criminals - events, A, config.grid.boundary);

    sigma = config.model.arrest_probability;
    generation_intensity = generation_ratio * generation_scale .* ...
        (1 - sigma) .* exp(-pi) * dt;
    switch config.model.generation_distribution
        case "bernoulli"
            generated = binornd(ones(ny, nx), -expm1(-generation_intensity));
        case "poisson"
            generated = poissrnd(generation_intensity);
    end

    diffusion = neighbor_sum(B, config.grid.boundary) - 4 * B;
    next_B = (B + 0.25 * config.model.eta * diffusion) * (1 - dt) ...
        + event_increment * events;
    if any(~isfinite(next_B), 'all') || min(next_B, [], 'all') < -1e-12
        error('urbancrime:AttractivenessPositivity', ...
            'Attractiveness update lost positivity.');
    end

    realized = event_unit * events;
    if config.policing.strategy == "delayed"
        busy_probability = -expm1(-sigma * criminals .* probability);
        staying = binornd(police, busy_probability);
        police_arrivals = move_population( ...
            police - staying, H, config.grid.boundary);
        next_police = staying + police_arrivals;
        next_H = H + dt / config.policing.tau * (realized - H);
    else
        next_police = police;
        next_H = zeros(ny, nx);
    end

    if sum(next_police, 'all') ~= sum(police, 'all')
        error('urbancrime:PoliceConservation', ...
            'Police movement failed to conserve the integer budget.');
    end
    if any(~isfinite(next_H), 'all') || min(next_H, [], 'all') < -1e-12
        error('urbancrime:InformationPositivity', ...
            'Delayed information update lost positivity.');
    end

    B = max(next_B, 0);
    criminals = criminal_arrivals + generated;
    police = next_police;
    H = max(next_H, 0);
    last_events = events;
    last_realized = realized;
    total_events = total_events + sum(events, 'all');
    total_generated = total_generated + sum(generated, 'all');

    slot = snapshot_slot(step + 1);
    if slot > 0
        save_snapshot(slot);
    end
    slot = metric_slot(step + 1);
    if slot > 0
        save_metrics(slot);
    end
end

result.summary.initial_criminals = initial_criminals;
result.summary.final_criminals = sum(criminals, 'all');
result.summary.total_generated = total_generated;
result.summary.total_events = total_events;
result.summary.initial_police_agents = initial_police;
result.summary.final_police_agents = sum(police, 'all');
result.config = config;

    function save_snapshot(slot)
        current_A = ast + B;
        current_rho = rho_unit * criminals;
        current_pi = pi_unit * police;
        result.A(:, :, slot) = current_A;
        result.rho(:, :, slot) = current_rho;
        result.events(:, :, slot) = last_events;
        result.pi(:, :, slot) = current_pi;
        result.H(:, :, slot) = H;
        result.expected_S(:, :, slot) = current_rho .* current_A .* exp(-current_pi);
        result.realized_S(:, :, slot) = last_realized;
        result.criminal_counts(:, :, slot) = criminals;
        result.police_counts(:, :, slot) = police;
    end

    function save_metrics(slot)
        current_A = ast + B;
        current_rho = rho_unit * criminals;
        current_pi = pi_unit * police;
        expected = current_rho .* current_A .* exp(-current_pi);
        result.metrics.mean_A(slot) = mean(current_A, 'all');
        result.metrics.mean_rho(slot) = mean(current_rho, 'all');
        result.metrics.mean_pi(slot) = mean(current_pi, 'all');
        result.metrics.mean_H(slot) = mean(H, 'all');
        result.metrics.mean_expected_S(slot) = mean(expected, 'all');
        result.metrics.mean_realized_S(slot) = mean(last_realized, 'all');
        result.metrics.cumulative_events(slot) = total_events;
    end
end

function values = expand_field(value, ny, nx, name)
validateattributes(value, {'numeric'}, {'real', 'finite', 'nonnegative'}, '', name);
if isscalar(value)
    values = repmat(double(value), ny, nx);
elseif isequal(size(value), [ny, nx])
    values = double(value);
else
    error('urbancrime:FieldShape', '%s must be scalar or grid-sized.', name);
end
end

function arrivals = move_population(movers, preference, boundary)
[ny, nx] = size(movers);
arrivals = zeros(ny, nx);
for i = 1:ny
    for j = 1:nx
        count = movers(i, j);
        if count == 0
            continue
        end
        neighbors = neighbor_indices(i, j, ny, nx, boundary);
        weights = zeros(1, 4);
        for direction = 1:4
            weights(direction) = preference( ...
                neighbors(1, direction), neighbors(2, direction));
        end
        if sum(weights) <= 0
            probabilities = 0.25 * ones(1, 4);
        else
            probabilities = weights / sum(weights);
        end
        draw = mnrnd(count, probabilities);
        for direction = 1:4
            ii = neighbors(1, direction);
            jj = neighbors(2, direction);
            arrivals(ii, jj) = arrivals(ii, jj) + draw(direction);
        end
    end
end
end

function total = neighbor_sum(values, boundary)
[ny, nx] = size(values);
total = zeros(ny, nx);
for i = 1:ny
    for j = 1:nx
        neighbors = neighbor_indices(i, j, ny, nx, boundary);
        for direction = 1:4
            total(i, j) = total(i, j) + values( ...
                neighbors(1, direction), neighbors(2, direction));
        end
    end
end
end

function indices = neighbor_indices(i, j, ny, nx, boundary)
% Order: west, south, east, north.
rows = [i, i - 1, i, i + 1];
cols = [j - 1, j, j + 1, j];
switch boundary
    case "periodic"
        rows(rows == 0) = ny;
        rows(rows > ny) = 1;
        cols(cols == 0) = nx;
        cols(cols > nx) = 1;
    case "no_flow"
        rows(rows == 0) = i + 1;
        rows(rows > ny) = i - 1;
        cols(cols == 0) = j + 1;
        cols(cols > nx) = j - 1;
end
indices = [rows; cols];
end

function validate_config(config)
required = {'grid', 'time', 'scaling', 'model', 'initial', 'policing', 'output'};
for index = 1:numel(required)
    if ~isfield(config, required{index})
        error('urbancrime:MissingConfig', 'Missing config.%s.', required{index});
    end
end
if config.grid.nx < 2 || config.grid.ny < 2 || config.grid.spacing <= 0
    error('urbancrime:InvalidGrid', 'Grid dimensions must be at least two.');
end
if ~ismember(config.grid.boundary, ["no_flow", "periodic"])
    error('urbancrime:InvalidBoundary', 'Boundary must be no_flow or periodic.');
end
raw_steps = (config.time.final - config.time.start) / config.time.step;
if config.time.step <= 0 || abs(raw_steps - round(raw_steps)) > 1e-8
    error('urbancrime:InvalidTime', 'Time interval must be divisible by time.step.');
end
if config.time.step > 1 || config.time.save_every < 1
    error('urbancrime:InvalidTime', 'Invalid time step or save interval.');
end
if config.model.eta < 0 || config.model.eta > 1
    error('urbancrime:InvalidEta', 'model.eta must lie in [0,1].');
end
if config.model.arrest_probability < 0 || config.model.arrest_probability > 1
    error('urbancrime:InvalidSigma', 'arrest_probability must lie in [0,1].');
end
if config.policing.strategy == "delayed"
    if config.policing.tau <= 0 || config.time.step > config.policing.tau
        error('urbancrime:InvalidDelay', 'Delayed policing requires tau >= time.step.');
    end
elseif config.policing.strategy ~= "none"
    error('urbancrime:InvalidStrategy', 'Strategy must be none or delayed.');
end
steps = (config.output.snapshot_times - config.time.start) / config.time.step;
if any(config.output.snapshot_times < config.time.start) || ...
        any(config.output.snapshot_times > config.time.final) || ...
        any(abs(steps - round(steps)) > 1e-8)
    error('urbancrime:InvalidSnapshots', 'Snapshot times must lie on time steps.');
end
end
