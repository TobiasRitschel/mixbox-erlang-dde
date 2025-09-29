%% Test numerical simulation with nonstiff method
% Clear command window
clc;

% Clear all variables
clear all;

% Close figures
close all;

% Remove added paths
restoredefaultpath;

%% Add paths
% Load library
run('../../../load_library');

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

% Set default interpreter
set(groot, 'DefaultTextInterpreter',            'Latex');
set(groot, 'DefaultAxesTickLabelInterpreter',   'Latex');
set(groot, 'DefaultLegendInterpreter',          'Latex');

%% System
% Random number generator seed
rng(0);

% Number of states and delayed quantities
nx = 2;
nz = 2;

% Parameters
kappa = 1;
K     = 1;

% Right-hand side function
f = @logistic_right_hand_side_function_multikernel;

% Delay function
h = @logistic_delay_function;

% Highest polynomial order in kernel
M = [9; 5];

% Forgetting factor parameter
a = [5; 2.5];

% Coefficients in kernel
c = abs(randn(2, max(M)+1));
c = c./sum(c, 2);

% Set "virtual" coefficients to zero
c(1:max(M)+1 > M+1) = 0;

% Number of time steps in simulation
N = 20000;

%% Simulate DDE transformed to ODEs using the LCT
% Parameter structs
p.kappa = kappa;
p.K     = K;

p.nx = nx;
p.nz = nz;

p.f = f;
p.h = h;

% Initial conditions
x0 = [0.9; 0.6];

% Time interval (simulation)
tspan = linspace(0, 20, N+1);

%% Simulate DDE using numerical approach for non-stiff systems
% Kernel
alpha = @(t) evaluateKernel(t, c, a);

% Memory horizon
Nh = N;

% Start timing
cpu_matlab_id = tic;

% Simulate system numerically
[T, X, Z, R] = numerical_simulation_nonstiff(f, tspan, x0, p, alpha, h, Nh);

% Stop timing
cpu_matlab = toc(cpu_matlab_id);

% Print results
fprintf('CPU time (Matlab): %.3f\n', cpu_matlab);

% Start timing
cpu_mex_id = tic;

% Simulate system numerically
[T_mex, X_mex, Z_mex, R_mex] = numerical_simulation_nonstiff_mex(f, tspan, x0, p, alpha, h, Nh);

% Stop timing
cpu_mex = toc(cpu_mex_id);

% Print results
fprintf('CPU time (MEX):    %.3f\n', cpu_mex);

%% Visualize simulations
% Create figure
figure(1);

% Select subplot
subplot(211);

% Visualize simulation obtained with Matlab routine
plot(T, X);

% Keep adding to the plot
hold on;

% Visualize simulation obtained with MEX routine
plot(T_mex, X_mex, '--');

% Stop adding to the plot
hold off;

% Axis limits
xlim(T([1, end]));

% Add legend
legend('x1 (Matlab)', 'x2 (Matlab)', 'x1 (MEX)', 'x2 (MEX)', 'Location', 'SouthEast');

% Select subplot
subplot(212);

% Visualize difference
plot(T, X - X_mex);

% Axis limits
xlim(T([1, end]));

% Add legend
legend('x1', 'x2', 'Location', 'SouthEast');