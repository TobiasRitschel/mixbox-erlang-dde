function [c, a, obj, exitflag, output] = identify_kernel(t, M, alpha, varargin)

% Options
opts = [];
if(nargin > 3), opts = varargin{1}; end

% fmincon options
fmincon_opts = optimoptions('fmincon');
if(~isempty(opts) && isfield(opts, 'fmincon_opts') && ~isempty(opts.fmincon_opts))
    fmincon_opts = opts.fmincon_opts;
end

% Use MEX routines?
UseMEXRoutines = true;
if(~isempty(opts) && isfield(opts, 'UseMEXRoutines'))
    UseMEXRoutines = opts.UseMEXRoutines;
end

% Use long or double in kernel computations?
UseLongKernelComputations = false;
if(~isempty(opts) && isfield(opts, 'UseLongKernelComputations'))
    UseLongKernelComputations = opts.UseLongKernelComputations;
end

switch UseMEXRoutines
    case true
        switch UseLongKernelComputations
            case true
                % MEX routines with long kernel computations
                objective_gradient_routine = @least_squares_kernel_objective_long_mex;
                hessian_routine            = @least_squares_kernel_objective_hessian_long_mex;
            case false
                % MEX routines with double kernel computations
                objective_gradient_routine = @least_squares_kernel_objective_mex;
                hessian_routine            = @least_squares_kernel_objective_hessian_mex;
        end
    case false
        % Matlab routines with double kernel computations
        objective_gradient_routine = @least_squares_kernel_objective;
        hessian_routine            = @least_squares_kernel_objective_hessian;
end

% Evaluate true kernel
alphameas = alpha(t, varargin{2:end});

% Overwrite specific fmincon settings
fmincon_opts = optimoptions(fmincon_opts, ...
    'SpecifyObjectiveGradient', true,          ...
    'HessianFcn',               hessian_routine);

%% Estimate kernel parameters
% Initial guess of interpolation coefficients and forgetting rate
c0 = ones(1, M+1)/(M+1);
a0 = 20;

% Collect initial guesses
theta0 = [c0, a0];

% Inequality constraint system matrix and right-hand side
A = [];
B = [];

% Equality constraint system matrix and right-hand side
Aeq = [ones(1, M+1), 0];
Beq = 1;

% Upper and lower bounds
ub = [ ones(M+1, 1); inf];
lb = [zeros(M+1, 1);   1];

% Nonlinear constraint function
nonlcon = [];

% Estimate the interpolation coefficients
[thetaest, obj, exitflag, output, lambda, grad, hessian] = ...
    fmincon(objective_gradient_routine, ...
    theta0, A, B, Aeq, Beq, lb, ub, nonlcon, fmincon_opts, ...
    t, alphameas, opts); %#ok

% Parameters in approximate (mixed Erlang) kernel
c = thetaest(1:end-1);
a = thetaest(end);