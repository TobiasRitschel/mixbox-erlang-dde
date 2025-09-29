%% Test Jacobian of right-hand side of distributed delay differential equation
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
% Add path to source code
run('../load_library');

% Load finite difference tool
run('../ext/FiniteDifferenceTool/LoadLibrary');

% Add path to utility functions
addpath(genpath(fullfile(pwd, './util')));

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

%% Load data
% Load data
load('./data/logistic_equation');

% Number of states
nx = numel(x0);

%% Finite difference test w. arbitrary parameters
% Plot results
PlottingOn = true;

% Finite differences
NoTests = 2e1;
eps = 10.^(linspace(-12, -2, NoTests));

% Error function
ErrorFunction = @(df, dfFD) abs(df(:) - dfFD(:));
% ErrorFunction = @(df, dfFD) abs(df(:) - dfFD(:))./(1 + abs(dfFD(:)));
% ErrorFunction = @(df, dfFD) abs(df(:) - dfFD(:))./(1e-8 + abs(dfFD(:)));
% ErrorFunction = @(df, dfFD) abs(df(:) - dfFD(:))./(norm(dfFD(:)) + abs(dfFD(:)));

% Number of parameters to be estimated (M+1 coefficients and 1 exponent)
p.ntheta = M+2 + nx + p.np;

% Estimation options
opts.simulator  = simulator;
opts.odeopts    = odeopts;

% Evaluate initial value of delayed quantities
z0 = h(x0, p);

% Initial values of all delayed quantities
Z0 = repmat(z0(:), M+1, 1);

% Initial condition
inp0 = [x0(:); Z0(:)];

% Add random perturbation
inp0 = inp0 + 1e2*randn(size(inp0));

% Carry out finite difference test
edf = FastFiniteDifferenceTest(@right_hand_side_ddde_wrapper, ...
    inp0, eps, PlottingOn, ErrorFunction, ...
    [], p);