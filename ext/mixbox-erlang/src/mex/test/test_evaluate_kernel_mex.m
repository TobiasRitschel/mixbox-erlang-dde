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

%% Run test
% Random number generator seed
rng(0);

% Number of coefficients
M = [150; 115];

% Number of kernels
nz = numel(M);

% Time
t = abs(randn(1, 1e3));

% Initial guess of interpolation coefficients and forgetting rate
c = abs(randn(nz, max(M)+1)/(max(M)+1));
a = [20; 10];

% Indices of "virtual" coefficients
idx = 1:max(M)+1 > M+1;

% Set "virtual" coefficients to zero
c(idx) = 0;

% Normalize
c = c./sum(c, 2);

% Evaluate kernel using Matlab routine
[alpha, dalpha, d2alpha] = evaluateKernel(t, c, a);

% Mex routines
evaluate_kernel_mex_routines = {    ...
    @evaluate_kernel_mex,           ...
    @evaluate_kernel_fast_mex,      ...
    @evaluate_kernel_faster_mex };

for idx = 1:numel(evaluate_kernel_mex_routines)
    % Select MEX routine
    evaluateKernelMex = evaluate_kernel_mex_routines{idx};

    % Print message to the user
    fprintf('\n\n --- %-s --- \n\n\n', func2str(evaluateKernelMex));

    %% Compare results with Matlab routine
    % Evaluate kernel using MEX routine
    [alpha_mex, dalpha_mex, d2alpha_mex] = evaluateKernelMex(t, c, a);

    % Display results
    fprintf('||alpha   - alpha_mex  || = %.3g\n', norm(alpha  (:) - alpha_mex  (:)));
    fprintf('||dalpha  - dalpha_mex || = %.3g\n', norm(dalpha (:) - dalpha_mex (:)));
    fprintf('||d2alpha - d2alpha_mex|| = %.3g\n', norm(d2alpha(:) - d2alpha_mex(:)));

    %% Compare computation time
    % Number of repetitions
    M = 1e2;

    for number_of_outputs = 1:3
        % Print message
        fprintf('\nNumber of outputs: %d\n', number_of_outputs);

        % Start timing
        cpu_matlab_id = tic;

        for i = 1:M
            % Evaluate kernel using Matlab routine
            [out{1:number_of_outputs}] = evaluateKernel(t, c, a);
        end

        % Stop timing
        cpu_matlab = toc(cpu_matlab_id);

        % Print results
        fprintf('CPU time (Matlab): %.3f\n', cpu_matlab);

        % Start timing
        cpu_mex_id = tic;

        for i = 1:M
            % Evaluate kernel using Matlab routine
            [out{1:number_of_outputs}] = evaluateKernelMex(t, c, a);
        end

        % Stop timing
        cpu_mex = toc(cpu_mex_id);

        % Print results
        fprintf('CPU time (MEX):    %.3f\n', cpu_mex);

        % Clear outputs
        clear out;
    end
end