%% Test Jacobian of right-hand side function of the logistic equation
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

%% System
% Parameters
kappa = 4;
K     = 1;

% Right-hand side function
f = @logistic_right_hand_side_function;

% Delay function
h = @logistic_delay_function;

% Parameter structs
p.kappa = kappa;
p.K     = K;

% Initial conditions
x0 = 0.9;
z0 = h(x0, p);

% Number of states and delayed quantities
p.nx = numel(x0);
p.nz = numel(z0);

%% Finite difference test w. arbitrary parameters
% Plot results
PlottingOn = true;

% Finite differences
NoTests = 2e1;
eps = 10.^(linspace(-6, 12, NoTests));

% Error function
ErrorFunction = @(df, dfFD) abs(df(:) - dfFD(:));
% ErrorFunction = @(df, dfFD) abs(df(:) - dfFD(:))./(1 + abs(dfFD(:)));
% ErrorFunction = @(df, dfFD) abs(df(:) - dfFD(:))./(1e-8 + abs(dfFD(:)));
% ErrorFunction = @(df, dfFD) abs(df(:) - dfFD(:))./(norm(dfFD(:)) + abs(dfFD(:)));

% Carry out finite difference test
edf = FastFiniteDifferenceTest(@right_hand_side_function_wrapper, ...
    [x0; z0], eps, PlottingOn, ErrorFunction, ...
    p, f);