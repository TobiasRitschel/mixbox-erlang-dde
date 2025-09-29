%% Simulate the logistic equation with a distributed delay
% Clear command window
clc;

% Clear all variables
clear all;

% Close figures
close all;

% Remove added paths
restoredefaultpath;

% Add path to source code
run('../load_library');

% Add external libraries
run('../ext/mixbox-erlang/load_library.m');

%% Formatting
% Save figures?
SaveFigures = true;

% Format of figures
Format = 'eps';

% Font size
fs = 18;

% Line width
lw = 3;

% Marker size
ms = 20;

% Aspect ratio
ar = [3, 1, 1];

% Figure position
fig_pos = [1, 2, 8, 6];

% Set default figure settings
set_default_figure_settings(fs, lw, ms, fig_pos);

%% System
% Parameters
kappa = 4;
K     = 1;

% Delay function
h = @logistic_delay_function;

% Simulator
simulator = @ode45;

% Initial and final time
t0 =  0;
tf = 24;

%% True solution and source term
% Time dilation parameter
gamma = 10;

% True solution
xtrue = @(t) 1 + exp(-(t/gamma).^2);

% Derivative of true solution
dxtrue = @(t) -2*t/gamma^2.*exp(-(t/gamma).^2);

% True delayed state
ztrue = @(t) 1 + gamma/sqrt(gamma^2 + 1)*exp(-t.^2/(gamma^2 + 1)).*(1 + erf(t/(gamma*sqrt(gamma^2 + 1))));

% Source term
p.Q = @(t, p) dxtrue(t) - p.kappa*xtrue(t).*(1 - ztrue(t)/p.K);

% Modified right-hand side function
f = @logistic_right_hand_side_function_with_source;

%% Simulation scenario
% Parameter structs
p.kappa = kappa;
p.K     = K;

% Initial conditions
x0 = xtrue(0);
z0 = h(x0, p);

% Number of time steps in numerical simulation (with Euler's implicit method)
Nsim = (tf - t0)*10^3 + 1;

% Time interval (simulation)
tspan = linspace(t0, tf, Nsim); % [mo]

% Number of states and delayed quantities
nx = numel(x0);
nz = numel(z0);

% Variable dimensions
p.nx = nx;
p.nz = nz;

% Function handles
p.f = f;
p.h = h;

% Tolerance
abstol = 1e-12;
reltol = 1e-12;

% Options
odeopts = odeset('AbsTol', abstol, 'RelTol', reltol, ...
    'Jacobian', @jacobian_ddde);

%% Kernel
% Store in regular array
alpha = @(t) 2/sqrt(pi)*exp(-t.^2);
beta  = @(t) erf(t);

%% Simulate using kernel approximation and the linear chain trick
% Tolerance
OptimalityTolerance = 1e-10;
StepTolerance       = 1e-20; % Just set it very low to disable it

% Maximum number of iterations
MaxIterations = 3e3;

% Bisection function tolerances
opts.bisection_tol = 1e-15;

% Optimizer options
opts.fmincon_opts = optimoptions('fmincon',          ...
    'Algorithm',                'interior-point',    ...
    'OptimalityTolerance',      OptimalityTolerance, ...
    'StepTolerance',            StepTolerance,       ...
    'MaxIterations',            MaxIterations,       ...
    'Display',                  'Iter');

% Tolerance for main support
epsilon = 1e-14;

% Number of measurements
N = 1000;

% Identify domain using bisection
th = identify_domain(beta, epsilon, opts);

% Measurement points
tmeas = linspace(0, th, N+1);

% Use long kernel computations
opts.UseLongKernelComputations = true;

% Kernel parameter
M = 2^7 - 1;

% Start timing
cpu_fmincon_id = tic;

% Identify kernel
[chat, ahat] = identify_kernel(tmeas, M, alpha, opts);

% Stop timing
cpu_time_fmincon = toc(cpu_fmincon_id);

% Display M and computation time
fprintf('M = %d\n', M);
fprintf('CPU time (s): %5.1f\n', cpu_time_fmincon);

% Report results
fprintf('\n');
fprintf('M%2s ', '');
fprintf('a%8s ', '');
fprintf('c%-6d ', 0:M);
fprintf('\n');
for i = 1:size(chat, 1)
    fprintf('%3d ', M);
    fprintf('%7.3f ', [ahat(i), chat(i, 1:M+1)]);
    fprintf('\n');
end

% Visualize estimated kernel
plot_kernel(alpha, chat, ahat, M, th, @evaluateKernel, 'Erlang mixture approximation: Estimated hyperparameters');

%% Theoretical approximation of kernel
% Evaluate theoretical approximation of the kernel
[cbar, abar] = identify_kernel_theoretical(beta, th, M);

% Visualize estimated kernel
plot_kernel(alpha, cbar, abar, M, th, @evaluateKernel, 'Erlang mixture approximation: Theoretical hyperparameters');

% Show figure
drawnow;

%% Simulate using the mixed Erlang distribution approximation and the linear chain trick
% Copy parameters
pbar = p;

% System matrices
[p   .A, p   .B, p   .C] = createSystemMatrices(chat, ahat, M);
[pbar.A, pbar.B, pbar.C] = createSystemMatrices(cbar, abar, M);

% Number of time points used to evaluate the kernel error
Nd = 1e6;

% Initial condition for augmented system
inp0             = zeros(nx + sum(M) + nz, 1);
inp0_theoretical = zeros(nx + sum(M) + nz, 1);

inp0            (1:nx) = x0;
inp0_theoretical(1:nx) = x0;
for j = 1:M+1
    % Local index (order)
    m = j-1;

    % Compute initial memory state
    inp0            (nx + (j-1)*nz + (1:nz)) = evaluate_initial_memory_state(xtrue, ahat, m, t0, th, Nd);
    inp0_theoretical(nx + (j-1)*nz + (1:nz)) = evaluate_initial_memory_state(xtrue, abar, m, t0, th, Nd);
end

% Message to user
fprintf('\nSIMULATING WITH TWO ERLANG MIXTURE APPROXIMATIONS ... ');

% Simulate using the mixed erlang distribution approximation and the linear
% chain trick
[tlct, xlct] = simulator(@right_hand_side_ddde, tspan, inp0, odeopts, ...
    p);

[~, xlct_theoretical] = simulator(@right_hand_side_ddde, tspan, inp0_theoretical, odeopts, ...
    pbar);

% Message to user
fprintf('DONE\n\n');

% Errors
e_lct                    = sum((xlct            (2:end, 1:nx)' - xtrue(tspan(2:end))).^2.*diff(tspan), 2);
e_lct_theoretical        = sum((xlct_theoretical(2:end, 1:nx)' - xtrue(tspan(2:end))).^2.*diff(tspan), 2);

e_lct_kernel             = evaluate_kernel_error(chat, ahat, th, alpha, Nd, @evaluateKernel);
e_lct_kernel_theoretical = evaluate_kernel_error(cbar, abar, th, alpha, Nd, @evaluateKernel);

%% Simulate using Euler's method and quadrature
% Time step size
dt = diff(tspan([1, 2]));

% Memory horizon
Nh = (find(tspan >= tf, 1, 'First')-1) + 1;

% fsolve options
fsolve_opts = optimoptions('fsolve',    ...
    'Display',                  'none', ...
    'FunctionTolerance',        1e-12,  ...
    'OptimalityTolerance',      1e-08,  ...
    'SpecifyObjectiveGradient', true,   ...
    'CheckGradient',            false);

% Initial state
thist  = -fliplr(linspace(0, (Nh-1)*dt, Nh));
x0_num = xtrue(thist);

% Message to user
fprintf('SIMULATING WITH TWO REFERENCE METHODS ... ');

% Simulate system with actual kernel
[t_expl, x_expl] = numerical_simulation_nonstiff(f, tspan, x0_num, p, alpha, h, Nh             );
[t_impl, x_impl] = numerical_simulation_stiff   (f, tspan, x0_num, p, alpha, h, Nh, fsolve_opts);

% Message to user
fprintf('DONE\n\n');

%% Visualize simulation
% Create figure
figure('NumberTitle', 'Off', 'Name', 'Comparison');

% Select subplot
subplot(211);

% Visualize simulation obtained with Euler's method and quadrature
plot(t_impl, xtrue(t_impl), 'k');

% Reset color ordering
set(gca, 'ColorOrderIndex', 1);

% Keep adding to plot
hold on;

% Visualize simulation obtained with Euler's method and quadrature
plot(t_impl, x_impl, '--');

% Visualize simulation obtained with approximate kernel
plot(tlct, xlct(:, 1:nx), '--');

% Stop adding to plot
hold off;

% Axis limits
xlim(tspan([1, end]));

% Add legend
legend;

% Legend
legend('Exact', 'Num.', 'LCT');

% Select subplot
subplot(212);

% Visualize absolute difference
semilogy(tlct, abs(xtrue(tlct) - x_impl(:, Nh:end)'));

% Keep adding to plot
hold on;

% Visualize absolute difference
plot(tlct, abs(xtrue(tlct) - xlct(:, 1:nx)));

% Stop adding to plot
hold off;

% Axis limits
xlim(tspan([1, end]));

% Legend
legend('Num.', 'LCT');

%% Functions
function plot_kernel(alpha, cest, aest, M, th, evaluate_kernel, varargin)
% Title
title = '';
if(nargin > 6), title = varargin{1}; end

% Number of rows and columns
NumRows = 3;
NumCols = 1;

% Axis tick labels for bar chart of estimated variables
XTickLabels = cell(1, M + 2);

for l = 1:M + 1
    XTickLabels{l} = "$c_{" + num2str(l-1) + "}$";
end
XTickLabels{max(M, M) + 2} = "$a$";

% Times
tkern = linspace(0, th, 1e6);

% Estimated kernel
alphahat = evaluate_kernel(tkern, cest, aest);

% Create figure
figure('NumberTitle', 'Off', 'Name', title);

% Select subplot
subplot(NumRows, NumCols, 1);

% Visualize kernel
plot(tkern, alpha(tkern), 'DisplayName', 'True');

% Add more plots
hold on;

% Visualize kernel
plot(tkern, alphahat, '--', 'DisplayName', 'Estimated');

% Stop adding plots
hold off;

% Axis limits
xlim([0, th]);

% Axis labels
xlabel('Time [s]');

% Create legend
legend;

% Select subplot
subplot(NumRows, NumCols, 2);

% Visualize kernel
plot(tkern, alphahat - alpha(tkern), 'DisplayName', 'Error');

% Axis limits
xlim([0, th]);

% Axis labels
xlabel('Time [s]');

% Select subplot
subplot(NumRows, NumCols, 3);

% Visualize estimated parameter estimates
b = bar(cest, 'Hist');

% Set colors
b(1).FaceColor = [0.8500, 0.3250, 0.0980];

% Axis limits
xlim([0, M+2]);
% ylim([0, 1]);

% Axis ticks
xticks(1:16:M+2);

% Axis tick labels
xticklabels(XTickLabels);

% Axis labels
xlabel('Parameters');

% Remove tick lines
h = gca;
h.XAxis.TickLength = [0, 0];
end

function [f, dfx, dfz] = logistic_right_hand_side_function_with_source(t, x, z, p)

% Source term function
Q = p.Q;

if(nargout > 1)
    % Evaluate right-hand side function and Jacobians
    [f, dfx, dfz] = logistic_right_hand_side_function(t, x, z, p);
else
    % Evaluate right-hand side function
    f = logistic_right_hand_side_function(t, x, z, p);
end

% Add source term
f = f + Q(t, p);
end

function e = evaluate_kernel_error(chat, ahat, th, alpha, N, evaluate_kernel)

% Time points
t = linspace(0, th, N+1);

% Evaluate kernels
alphameas = alpha(t(1:end-1));
alphahat  = evaluate_kernel(t(1:end-1), chat, ahat);

% Approximate L1 error using a left rectangle rule
e = sum((alphameas - alphahat).^2.*diff(t), 2);
end

function zm0 = evaluate_initial_memory_state(x, a, m, t0, th, N)

% Time points
t = linspace(-th, 0, N);

% Copy time points
s = t;

% Evaluate subkernel
alpham = a^(m+1)/factorial(m)*(t0 - s).^m.*exp(-a*(t0 - s));

% Approximate integral using a right rectangle rule
zm0 = sum(x(t(2:end)).*alpham(2:end).*diff(t), 2);
end