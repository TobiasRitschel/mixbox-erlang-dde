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
mex ./src/mex/src/evaluate_kernel_mex.c;
mex ./src/mex/src/evaluate_kernel_fast_mex.c;
mex ./src/mex/src/evaluate_kernel_faster_mex.c;
mex ./src/mex/src/evaluate_kernel_faster_long_mex.c;
mex ./src/mex/src/least_squares_kernel_objective_mex.c;
mex ./src/mex/src/least_squares_kernel_objective_long_mex.c;
mex ./src/mex/src/least_squares_kernel_objective_hessian_mex.c;
mex ./src/mex/src/least_squares_kernel_objective_hessian_long_mex.c;

% Move files
movefile('./evaluate_kernel_mex.*',                             './src/mex/bin/', 'f');
movefile('./evaluate_kernel_fast_mex.*',                        './src/mex/bin/', 'f');
movefile('./evaluate_kernel_faster_mex.*',                      './src/mex/bin/', 'f');
movefile('./evaluate_kernel_faster_long_mex.*',                 './src/mex/bin/', 'f');
movefile('./least_squares_kernel_objective_mex.*',              './src/mex/bin/', 'f');
movefile('./least_squares_kernel_objective_long_mex.*',         './src/mex/bin/', 'f');
movefile('./least_squares_kernel_objective_hessian_mex.*',      './src/mex/bin/', 'f');
movefile('./least_squares_kernel_objective_hessian_long_mex.*', './src/mex/bin/', 'f');