function tests = test_matlab_abm
%TEST_MATLAB_ABM Fast checks for the maintained MATLAB implementation.
tests = functiontests(localfunctions);
end

function setupOnce(test_case)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root));
test_case.TestData.root = root;
end

function testPaperCaseMetadata(test_case)
m3as = matlab_abm_config("m3as", 3);
verifyEqual(test_case, m3as.model.eta, 0.03);
verifyEqual(test_case, m3as.output.snapshot_times, [0, 1.6, 10, 20, 200]);
extension = matlab_abm_config("police", 8);
verifyEqual(test_case, extension.policing.tau, 50.0);
verifyEqual(test_case, extension.output.snapshot_times, ...
    [0, 765, 768, 771, 774, 777, 780]);
end

function testPopulationRoundingAboveOneAgentPerCell(test_case)
rng(11, 'twister');
expected = 2.4 * ones(4, 3);
population = matlab_abm_integer_population(expected);
verifyEqual(test_case, sum(population, 'all'), ceil(sum(expected, 'all')));
verifyGreaterThanOrEqual(test_case, min(population, [], 'all'), 2);
verifyEqual(test_case, population, round(population));
end

function testNoPoliceSmoke(test_case)
config = small_config(matlab_abm_config("m3as", 1));
config.policing.strategy = "none";
config.initial.police_density = 0.0;
config.initial.information = 0.0;
result = matlab_abm_simulate(config);
verifyEqual(test_case, result.time, [0, 0.04, 0.10]);
verifyEqual(test_case, result.pi, zeros(size(result.pi)));
verifyEqual(test_case, result.H, zeros(size(result.H)));
verifyGreaterThanOrEqual(test_case, min(result.criminal_counts, [], 'all'), 0);
verifyEqual(test_case, result.criminal_counts, round(result.criminal_counts));
end

function testDelayedPoliceBudgetIsConserved(test_case)
config = small_config(matlab_abm_config("police", 1));
config.model.arrest_probability = 0.2;
result = matlab_abm_simulate(config);
totals = squeeze(sum(result.police_counts, [1, 2]));
verifyEqual(test_case, totals, repmat(totals(1), size(totals)));
verifyGreaterThanOrEqual(test_case, min(result.H, [], 'all'), 0);
verifyGreaterThanOrEqual(test_case, min(result.pi, [], 'all'), 0);
end

function config = small_config(config)
config.grid.nx = 6;
config.grid.ny = 5;
config.grid.spacing = 0.5;
config.time.start = 0.0;
config.time.final = 0.10;
config.time.step = 0.02;
config.time.save_every = 2;
config.output.snapshot_times = [0, 0.04, 0.10];
end
