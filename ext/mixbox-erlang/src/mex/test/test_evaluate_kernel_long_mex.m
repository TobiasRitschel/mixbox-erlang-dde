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
M = [500; 115];

% Number of kernels
nz = numel(M);

% Time
t = abs(randn(1, 1e3));

% Coefficients and rate parameter
c = abs(randn(nz, max(M)+1)/(max(M)+1));
a = [20; 10];

% Indices of "virtual" coefficients
idx = 1:max(M)+1 > M+1;

% Set "virtual" coefficients to zero
c(idx) = 0;

% Normalize
c = c./sum(c, 2);

% Evaluate kernel using Matlab routine
[alpha_double, dalpha_double, d2alpha_double] = evaluate_kernel_faster_mex(t, c, a);

% Evaluate kernel using MEX routine
[alpha_long, dalpha_long, d2alpha_long] = evaluate_kernel_faster_long_mex(t, c, a);

% Display results
fprintf('||alpha   (double) - alpha   (long)|| = %.3g\n', norm(alpha_double  (:) - alpha_long  (:)));
fprintf('||dalpha  (double) - dalpha  (long)|| = %.3g\n', norm(dalpha_double (:) - dalpha_long (:)));
fprintf('||d2alpha (double) - d2alpha (long)|| = %.3g\n', norm(d2alpha_double(:) - d2alpha_long(:)));

%% Compare computation time
% Number of repetitions
K = 1e2;

for number_of_outputs = 1:3
    % Print message
    fprintf('\nNumber of outputs: %d\n', number_of_outputs);

    % Start timing
    cpu_double_id = tic;

    for i = 1:K
        % Evaluate kernel using Matlab routine
        [out{1:number_of_outputs}] = evaluate_kernel_faster_mex(t, c, a);
    end

    % Stop timing
    cpu_double = toc(cpu_double_id);

    % Print results
    fprintf('CPU time (double): %.3f\n', cpu_double);

    % Start timing
    cpu_long_id = tic;

    for i = 1:K
        % Evaluate kernel using Matlab routine
        [out{1:number_of_outputs}] = evaluate_kernel_faster_long_mex(t, c, a);
    end

    % Stop timing
    cpu_long = toc(cpu_long_id);

    % Print results
    fprintf('CPU time (long):   %.3f\n', cpu_long);
end

%% Visualize difference
% Update coefficients
c = zeros(size(c));
c(1, M(1)+1) = 1;
c(2, M(2)+1) = 1;
c = abs(c)./sum(abs(c), 2);

% Update rate parameters
a = [100; 200];

% Times
t = linspace(0, 10, 1e3);

% Evaluate kernel using Matlab routine
alpha_double = evaluate_kernel_faster_mex(t, c, a);

% Evaluate kernel using MEX routine
alpha_long = evaluate_kernel_faster_long_mex(t, c, a);

% Create figure
figure(1);

% Select subplot
subplot(211);

% Plot kernels
plot(t, alpha_double, t, alpha_long, '--');

% Select subplot
subplot(212);

% Plot coefficients
plot(c');