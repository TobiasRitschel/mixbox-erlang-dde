function [c, a] = identify_kernel_theoretical(beta, th, M)

% Time step
dt = th/(M+1);

% Create time grid (used to evaluate coefficients)
tbar = 0:dt:th;

% Evaluate coefficients
c = diff(beta(tbar));

% Evaluate rate parameter
a = 1/dt;