%% Test MEX routine for evaluating a mixed Erlang kernel and its derivatives
% Clear command window
clc;

% Clear all variables
clear all;

% Close figures
close all;

% Remove added paths
restoredefaultpath;

% Reset graphical settings
reset(groot);

%% Add paths
% Load library
run('../../../load_library');

% Load finite difference tool
run('../../../ext/FiniteDifferenceTool/LoadLibrary');

%% Formatting
% Font size
fs = 12;

% Line width
lw = 3;

% Marker size
ms = 20;

% Set default font size
set(groot, 'DefaultAxesFontSize',   fs);

% Set default line widths
set(groot, 'DefaultLineLineWidth',  lw);
set(groot, 'DefaultStairLineWidth', lw);
set(groot, 'DefaultStemLineWidth',  lw);

% Set default marker size
set(groot, 'DefaultLineMarkerSize', ms);

%% Compare results with Matlab routine
% Random number generator seed
rng(0);

% Number of coefficients
M = 150;

% Number of kernels
nz = numel(M);

% Time
t = sort(abs(randn(1, 1e3)));

% Initial guess of interpolation coefficients and forgetting rate
c = abs(randn(nz, max(M)+1)/(max(M)+1));
a = 20;

% Measurements of the true kernel
alphameas = abs(randn(size(t)));
alphameas = alphameas/(sum(alphameas(1:end-1).*diff(t)));

% Normalize
c = c./sum(c, 2);

% Collect inputs
theta = [c, a];

% Evaluate kernel using Matlab routine
[phi, grad] = least_squares_kernel_objective        (theta,     t, alphameas, []);
hess        = least_squares_kernel_objective_hessian(theta, [], t, alphameas, []);

% Evaluate kernel using MEX routine
[phi_mex, grad_mex] = least_squares_kernel_objective_mex        (theta,     t, alphameas);
hess_mex            = least_squares_kernel_objective_hessian_mex(theta, [], t, alphameas);

% Display results
fprintf('||phi   - phi_mex || = %.3g\n', norm(phi (:) - phi_mex (:)));
fprintf('||grad  - grad_mex|| = %.3g\n', norm(grad(:) - grad_mex(:)));
fprintf('||hess  - hess_mex|| = %.3g\n', norm(hess(:) - hess_mex(:)));

%% Compare computation time
% Number of repetitions
M = 1e1;

for number_of_outputs = 1:2
    % Print message
    fprintf('\nNumber of outputs: %d\n', number_of_outputs);

    % Start timing
    cpu_matlab_id = tic;

    for i = 1:M
        % Evaluate kernel using Matlab routine
        [out{1:number_of_outputs}] = least_squares_kernel_objective(theta, t, alphameas, []);
    end

    % Stop timing
    cpu_matlab = toc(cpu_matlab_id);

    % Print results
    fprintf('CPU time (Matlab): %.3f\n', cpu_matlab);

    % Start timing
    cpu_mex_id = tic;

    for i = 1:M
        % Evaluate kernel using Matlab routine
        [out{1:number_of_outputs}] = least_squares_kernel_objective_mex(theta, t, alphameas);
    end

    % Stop timing
    cpu_mex = toc(cpu_mex_id);

    % Print results
    fprintf('CPU time (MEX):    %.3f\n', cpu_mex);

    % Clear outputs
    clear out;
end

% Print message
fprintf('\nNumber of outputs: 3\n');

% Start timing
cpu_matlab_id = tic;

for i = 1:M
    % Evaluate kernel using Matlab routine
    hess = least_squares_kernel_objective_hessian(theta, [], t, alphameas, []);
end

% Stop timing
cpu_matlab = toc(cpu_matlab_id);

% Print results
fprintf('CPU time (Matlab): %.3f\n', cpu_matlab);

% Start timing
cpu_mex_id = tic;

for i = 1:M
    % Evaluate kernel using Matlab routine
    hess_mex = least_squares_kernel_objective_hessian_mex(theta, [], t, alphameas);
end

% Stop timing
cpu_mex = toc(cpu_mex_id);

% Print results
fprintf('CPU time (MEX):    %.3f\n', cpu_mex);