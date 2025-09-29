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

%% Initialize
% Random number generator seed
rng(0);

% Orders
M = randi(10, 5, 1);

% Rate parameters
a = abs(randn(size(M, 1), 1));

% Coefficients
c = abs(randn(size(M, 1), max(M)+1));

% Indices of "virtual" coefficients
idx = 1:max(M)+1 > M+1;

% Set "virtual" coefficients to zero
c(idx) = 0;

% Normalize coefficients
c = c./sum(c, 2);

%% Compare sparse and dense versions
% Evaluate dense matrices
[Ad, Bd, Cd] = createSystemMatrices(c, a, M);

% Evaluate sparse matrices
[As, Bs, Cs] = createSystemMatrices(c, a, M, 'sparse');

% Report discrepancies
fprintf('|Ad - As| = %.2g\n', norm(Ad - full(As)));
fprintf('|Bd - Bs| = %.2g\n', norm(Bd - full(Bs)));
fprintf('|Cd - Cs| = %.2g\n', norm(Cd - full(Cs)));

% Create figure
figure(1);

% Clear figure
clf;

% Select subplot
subplot(321);

% Sparsity pattern of dense A
spy(Ad);

% Select subplot
subplot(322);

% Sparsity pattern of sparse A
spy(As);

% Select subplot
subplot(323);

% Sparsity pattern of dense B
spy(Bd);

% Select subplot
subplot(324);

% Sparsity pattern of sparse B
spy(Bs);

% Select subplot
subplot(325);

% Sparsity pattern of dense C
spy(Cd);

% Select subplot
subplot(326);

% Sparsity pattern of sparse C
spy(Cs);