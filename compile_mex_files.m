%% Compile the MEX functions located in ./src and move them to ./bin
% The compilation may complain that you are using a too new version of GCC,
% but this script intentionally does not specify a compiler to ensure that
% it runs on many different computers.

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

%% Load library
% Compile matricization programs
mex ./src/mex/src/numerical_simulation_nonstiff_mex.c;

% Move files
movefile('./numerical_simulation_nonstiff_mex.*', './src/mex/bin/', 'f');