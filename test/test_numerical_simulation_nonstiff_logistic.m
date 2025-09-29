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
% Add path to source code
run('../load_library');

% Load finite difference tool
run('../ext/FiniteDifferenceTool/LoadLibrary');

% Add external libraries
run('../ext/mixbox-erlang/load_library.m');

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
nx = 1;
nz = 1;

% Parameters
kappa = 1;
K     = 1;

% Right-hand side function
f = @logistic_right_hand_side_function;

% Delay function
h = @logistic_delay_function;

% Simulator
simulator = @ode45;

% Highest polynomial order in kernel
M = 9;

% Forgetting factor parameter
a = 5;

% Coefficients in kernel
c = abs(randn(1, M+1));
c = c./sum(c, 2);

% Number of time steps in simulation
N = 10000;

%% Simulate DDE transformed to ODEs using the LCT
% Parameter structs
p.kappa = kappa;
p.K     = K;

p.nx = nx;
p.nz = nz;

p.f = f;
p.h = h;

% System matrices
[p.A, p.B, p.C] = createSystemMatrices(c, a, M);

% Initial conditions
x0 = 0.9;
z0 = h(x0, p);

% Initial condition for augmented system
inp0 = zeros(M + 2*nz, 1);
inp0(1:nx      ) = x0;
inp0(  nx+1:end) = repelem(z0, M+1);

% Options
odeopts = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);

% Time interval (simulation)
tspan = linspace(0, 20, N+1);

% Simulate systems with approximate and actual kernels
[tlct, xlct] = simulator(@right_hand_side_ddde, tspan, inp0, odeopts, ...
    p);

%% Simulate DDE using numerical approach for non-stiff systems
% Kernel
alpha = @(t) evaluateKernel(t, c, a);

% Memory horizon
Nh = N;

% Simulate system numerically
[T, X, Z, R] = numerical_simulation_nonstiff(f, tspan, x0, p, alpha, h, Nh);

%% Visualize simulations
% Create figure
figure(1);

% Select subplot
subplot(211);

% Visualize simulation obtained with an ODE simulation and the LCT
plot(tlct, xlct(:, 1:nx));

% Keep adding to the plot
hold on;

% Visualize simulation obtained with numerical method (for non-stiff systems)
plot(T, X, '--');

% Stop adding to the plot
hold off;

% Axis limits
xlim(tspan([1, end]));

% Add legend
legend('x (LCT)', 'x (numerical)', 'Location', 'SouthEast');

% Select subplot
subplot(212);

% Visualize difference
plot(tlct, xlct(:, 1:nx) - X(:, Nh + (0:N))');

% Add legend
legend('Error', 'Location', 'SouthEast');